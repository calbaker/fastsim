"""
Calibration script for 2021_Hyundai_Sonata_Hybrid_Blue    
"""
from pathlib import Path

import numpy as np  # noqa: F401
import matplotlib.pyplot as plt  # noqa: F401
import seaborn as sns
import pandas as pd  # noqa: F401
import polars as pl  # noqa: F401
from typing import List, Dict
from pymoo.core.problem import StarmapParallelization
from copy import deepcopy

import fastsim as fsim
from fastsim import pymoo_api

mps_per_mph = 0.447
celsius_to_kelvin_offset = 273.15

# Initialize seaborn plot configuration
sns.set()

veh = fsim.Vehicle.from_file(Path(__file__).parent / "f3-vehicles/2021_Hyundai_Sonata_Hybrid_Blue.yaml")
veh_dict = veh.to_pydict()

sim_params_dict = fsim.SimParams.default().to_pydict()
sim_params_dict["trace_miss_opts"] = "AllowChecked"
sim_params = fsim.SimParams.from_pydict(sim_params_dict, skip_init=False)

# Obtain the data from
# https://nrel.sharepoint.com/:f:/r/sites/EEMSCoreModelingandDecisionSupport2022-2024/Shared%20Documents/FASTSim/DynoTestData?csf=1&web=1&e=F4FEBp
# and then copy it to the local folder below
cyc_folder_path = Path(__file__).parent / "dyno_test_data/2021 Hyundai Sonata Hybrid/Extended Datasets"
assert cyc_folder_path.exists()

# See 2021_Hyundai_Sonata_Hybrid_TestSummary_2022-03-01_D3.xlsx for cycle-level data
cyc_files: List[str] = [
    # The hot and cold cycles must have HVAC active!
    # - wide range of initial and ambient temperatures
    # - good signal quality -- somewhat subjective

    # HWY x2, hot (M155), HVAC active (B155)
    "62202004 Test Data.txt", 

    # US06 x2, hot, HVAC active
    "62202005 Test Data.txt",

    # UDDS x1, room temperature ambient
    "62201013 Test Data.txt",

    # HWY x2, room temperature ambient
    "62201014 Test Data.txt",

    # UDDSx2, 4 bag (FTP), cold start, in COLD (20°F) test cell, HVAC-AUTO-72°F, ECO drive mode
    "62202013 Test Data.txt",

    # UDDS, 2 bag, warm start, in COLD (20°F) test cell, HVAC-AUTO-72°F, ECO drive mode
    "62202014 Test Data.txt",

    # US06x2, 4 (split) bag, warm start, in COLD (20°F) test cell, HVAC-AUTO-72°F, ECO drive mode
    "62202016 Test Data.txt",

    # TODO: check for seat heater usage in cold cycles and account for that in model!
]
assert len(cyc_files) > 0
cyc_files: List[Path] = [cyc_folder_path / cyc_file for cyc_file in cyc_files]
print("\ncyc_files:\n", '\n'.join([cf.name for cf in cyc_files]), sep='')

# use random or manual selection to retain ~70% of cycles for calibration,
# and reserve the remaining for validation
cyc_files_for_cal: List[str] = [
    "62202004 Test Data.txt", 
    # "62202005 Test Data.txt",
    "62201013 Test Data.txt",
    "62201014 Test Data.txt",
    "62202013 Test Data.txt",
    # "62202014 Test Data.txt", 
    "62202016 Test Data.txt",
]
cyc_files_for_cal: List[Path] = [cyc_file for cyc_file in cyc_files if cyc_file.name in cyc_files_for_cal]
assert len(cyc_files_for_cal) > 0
print("\ncyc_files_for_cal:\n", '\n'.join([cf.name for cf in cyc_files_for_cal]), sep='')

time_column = "Time[s]_RawFacilities"
speed_column = "Dyno_Spd[mph]"

def df_to_cyc(df: pd.DataFrame) -> fsim.Cycle:
    cyc_dict = {
        "time_seconds": df[time_column].to_list(),
        "speed_meters_per_second": (df[speed_column] * mps_per_mph).to_list(),
        "temp_amb_air_kelvin": (df["Cell_Temp[C]"] + celsius_to_kelvin_offset).to_list(),
        # TODO: pipe solar load from `Cycle` into cabin thermal model
        # TODO: use something (e.g. regex) to determine solar load
        # see column J comments in 2021_Hyundai_Sonata_Hybrid_TestSummary_2022-03-01_D3.xlsx
        # "pwr_solar_load_watts": df[],
    }
    return fsim.Cycle.from_pydict(cyc_dict, skip_init=False)

def veh_init(cyc_file_stem: str, dfs: Dict[str, pd.DataFrame]) -> fsim.Vehicle:
    vd = deepcopy(veh_dict)
    # initialize SOC
    vd['pt_type']['HybridElectricVehicle']['res']['state']['soc'] = \
        dfs[cyc_file_stem]["HVBatt_SOC_high_precision_PCAN__per"].iloc[1] / 100
    assert 0 < vd['pt_type']['HybridElectricVehicle']['res']['state']['soc'] < 1, "\ninit soc: {}\nhead: {}".format(
        vd['pt_type']['HybridElectricVehicle']['res']['state']['soc'], dfs[cyc_file_stem]["HVBatt_SOC_high_precision_PCAN__per"].head())
    # initialize cabin temp
    vd['cabin']['LumpedCabin']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["Cabin_Temp[C]"][0] + celsius_to_kelvin_offset
    # initialize battery temperature to match cabin temperature because battery
    # temperature is not available in test data
    # Also, battery temperature has no effect in the HEV because efficiency data
    # does not go below 23*C and there is no active thermal management
    vd['pt_type']['HybridElectricVehicle']['res']['thrml']['RESLumpedThermal']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["Cabin_Temp[C]"][0] + celsius_to_kelvin_offset
    # initialize engine temperature
    vd['pt_type']['HybridElectricVehicle']['fc']['thrml']['FuelConverterThermal']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["engine_coolant_temp_PCAN__C"][0] + celsius_to_kelvin_offset
    return fsim.Vehicle.from_pydict(vd, skip_init=False)


dfs_for_cal: Dict[str, pd.DataFrame] = {
    # `delimiter="\t"` should work for tab separated variables
    cyc_file.stem: pd.read_csv(cyc_file, delimiter="\t") for cyc_file in cyc_files_for_cal
}
for key, df_for_cal in dfs_for_cal.items():
    # filter out "before" time
    df_for_cal = df_for_cal[df_for_cal[time_column] >= 0.0]
    # TODO: figure out if we should use an integrator for resampling rate vars
    # df_for_cal = df_for_cal.set_index(time_column)
    # df_for_cal = df_for_cal.resample("1s", origin="start").bfill()
    df_for_cal = df_for_cal[::10]
    df_for_cal.reset_index(inplace=True)
    dfs_for_cal[key] = df_for_cal
    
cycs_for_cal: Dict[str, fsim.Cycle] = {}
# populate `cycs_for_cal`
for (cyc_file_stem, df) in dfs_for_cal.items():
    cyc_file_stem: str
    df: pd.DataFrame
    cyc_dict_raw = df.to_dict()
    cyc_file_stem: str
    df: pd.DataFrame
    cycs_for_cal[cyc_file_stem] = df_to_cyc(df)

sds_for_cal: Dict[str, fsim.SimDrive] = {}
# populate `sds_for_cal`
for (cyc_file_stem, cyc) in cycs_for_cal.items():
    cyc_file_stem: str
    cyc: fsim.Cycle
    # NOTE: maybe change `save_interval` to 5
    veh = veh_init(cyc_file_stem, dfs_for_cal)
    sds_for_cal[cyc_file_stem] = fsim.SimDrive(veh, cyc, sim_params).to_pydict()

cyc_files_for_val: List[Path] = list(set(cyc_files) - set(cyc_files_for_cal))
assert len(cyc_files_for_val) > 0
print("\ncyc_files_for_val:\n", '\n'.join([cf.name for cf in cyc_files_for_val]), sep='')

dfs_for_val: Dict[str, pd.DataFrame] = {
    # `delimiter="\t"` should work for tab separated variables
    cyc_file.stem: pd.read_csv(cyc_file, delimiter="\t") for cyc_file in cyc_files_for_val
}
for key, df_for_val in dfs_for_val.items():
    # filter out "before" time
    df_for_val = df_for_val[df_for_val[time_column] >= 0.0]
    # TODO: figure out if we should use an integrator for resampling rate vars
    # df_for_val = df_for_val.set_index(time_column)
    # df_for_val = df_for_val.resample("1s", origin="start").bfill()
    df_for_val = df_for_val[::10]
    df_for_val.reset_index(inplace=True)
    dfs_for_val[key] = df_for_val

cycs_for_val: Dict[str, fsim.Cycle] = {}
# populate `cycs_for_val`
for (cyc_file_stem, df) in dfs_for_val.items():
    cyc_file_stem: str
    df: pd.DataFrame
    cycs_for_val[cyc_file_stem] = df_to_cyc(df)

sds_for_val: Dict[str, fsim.SimDrive] = {}
# populate `sds_for_val`
for (cyc_file_stem, cyc) in cycs_for_val.items():
    cyc_file_stem: str
    cyc: fsim.Cycle
    veh = veh_init(cyc_file_stem, dfs_for_val)
    sds_for_val[cyc_file_stem] = fsim.SimDrive(veh, cyc, sim_params).to_pydict()

# Setup model objectives
## Parameter Functions
def new_em_eff_max(sd_dict, new_eff_max) -> Dict:
    """
    Set `new_eff_max` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['HybridElectricVehicle']['em'])
    em.__eff_fwd_max = new_eff_max
    sd_dict['veh']['pt_type']['HybridElectricVehicle']['em'] = em.to_pydict()
    return sd_dict

def new_em_eff_range(sd_dict, new_eff_range) -> Dict:
    """
    Set `new_eff_range` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['HybridElectricVehicle']['em'])
    em.__eff_fwd_range = new_eff_range
    sd_dict['veh']['pt_type']['HybridElectricVehicle']['em'] = em.to_pydict()
    return sd_dict

def new_fc_eff_max(sd_dict, new_eff_max) -> Dict:
    """
    Set `new_eff_max` in `FuelConverter`
    """
    fc = fsim.FuelConverter.from_pydict(sd_dict['veh']['pt_type']['HybridElectricVehicle']['fc'])
    fc.__eff_max = new_eff_max
    sd_dict['veh']['pt_type']['HybridElectricVehicle']['fc'] = fc.to_pydict()
    return sd_dict

def new_fc_eff_range(sd_dict, new_eff_range) -> Dict:
    """
    Set `new_eff_range` in `FuelConverter`
    """
    fc = fsim.FuelConverter.from_pydict(sd_dict['veh']['pt_type']['HybridElectricVehicle']['fc'])
    fc.__eff_range = new_eff_range
    sd_dict['veh']['pt_type']['HybridElectricVehicle']['fc'] = fc.to_pydict()
    return sd_dict

def get_mod_soc(sd_dict):
    return np.array(sd_dict['veh']['pt_type']['HybridElectricVehicle']['res']['history']['soc'])

def get_exp_soc(df):
    return df['HVBatt_SOC_high_precision_PCAN__per'] / 100

save_path = Path(__file__).parent / "pymoo_res" / Path(__file__).stem
save_path.mkdir(exist_ok=True, parents=True)

## Model Objectives
cal_mod_obj = pymoo_api.ModelObjectives(
    models = sds_for_cal,
    dfs = dfs_for_cal,
    obj_fns=(
        (
            get_mod_soc,
            get_exp_soc
        ),
        # TODO: add objectives for:
        # - achieved and cycle speed
        # - engine fuel usage 
        # - battery temperature -- BEV only, if available
        # - engine temperature
        # - cabin temperature
        # - HVAC power for cabin, if available
    ),
    param_fns=(
        new_em_eff_max,
        new_em_eff_range,
        new_fc_eff_max,
        # new_fc_eff_range, 
        # TODO: make sure this has functions for modifying
        # - cabin thermal
        #     - thermal mass
        #     - length
        #     - htc to amb when stopped
        #     - set width from vehicle specs -- no need to calibrate
        # - battery thermal -- not necessary for HEV because battery temperature has no real effect
        #     - thermal mass
        #     - convection to ambient
        #     - convection to cabin
        # ## HEV specific stuff
        # - HVAC PID controls for cabin (not for battery because Sonata has
        #   passive thermal management, but make sure to do battery thermal
        #   controls for BEV)
        # - engine thermal
        #     - thermal mass
        #     - convection to ambient when stopped
        #     - diameter
    ),
    # must match order and length of `params_fns`
    bounds=(
        (0.80, 0.99),
        (0.1, 0.6),
        (0.32, 0.45),
        # (0.0, 0.45),
    ),
    verbose=False,    
)

val_mod_obj = pymoo_api.ModelObjectives(
    models = sds_for_val,
    dfs = dfs_for_val,
    obj_fns=(
        (
            get_mod_soc,
            get_exp_soc
        ),
        # TODO: add objectives for:
        # - achieved and cycle speed
        # - engine fuel usage 
        # - battery temperature -- BEV only, if available
        # - engine temperature
        # - cabin temperature
        # - HVAC power for cabin, if available
    ),
    param_fns=(
        new_em_eff_max,
        new_em_eff_range,
        new_fc_eff_max,
        # new_fc_eff_range, 
        # TODO: make sure this has functions for modifying
        # - cabin thermal
        #     - thermal mass
        #     - length
        #     - htc to amb when stopped
        #     - set width from vehicle specs -- no need to valibrate
        # - battery thermal -- not necessary for HEV because battery temperature has no real effect
        #     - thermal mass
        #     - convection to ambient
        #     - convection to cabin
        # ## HEV specific stuff
        # - HVAC PID controls for cabin (not for battery because Sonata has
        #   passive thermal management, but make sure to do battery thermal
        #   controls for BEV)
        # - engine thermal
        #     - thermal mass
        #     - convection to ambient when stopped
        #     - diameter
    ),
    # must match order and length of `params_fns`
    bounds=(
        (0.80, 0.99),
        (0.1, 0.6),
        (0.32, 0.45),
        # (0.0, 0.45),
    ),
    verbose=False,    
)

em_eff_fwd_max = fsim.ElectricMachine.from_pydict(veh_dict['pt_type']['HybridElectricVehicle']['em'], skip_init=False).eff_fwd_max 
em_eff_fwd_range = fsim.ElectricMachine.from_pydict(veh_dict['pt_type']['HybridElectricVehicle']['em'], skip_init=False).eff_fwd_range 
fc_eff_max = fsim.FuelConverter.from_pydict(veh_dict['pt_type']['HybridElectricVehicle']['fc'], skip_init=False).eff_max 
# print("Verifying that model responds to input parameter changes by individually perturbing parameters")
# baseline_errors = cal_mod_obj.get_errors(
#     cal_mod_obj.update_params([
#         em_eff_fwd_max,
#         em_eff_fwd_range,
#         fc_eff_max,
#         # veh_dict['pt_type']['HybridElectricVehicle']['fc'],
#     ])
# )
# param0_perturb = cal_mod_obj.get_errors(
#     cal_mod_obj.update_params([
#         em_eff_fwd_max + 0.05,
#         em_eff_fwd_range,
#         fc_eff_max,
#         # veh_dict['pt_type']['HybridElectricVehicle']['fc'],
#     ])
# )
# assert list(param0_perturb.values()) != list(baseline_errors.values())
# param1_perturb = cal_mod_obj.get_errors(
#     cal_mod_obj.update_params([
#         em_eff_fwd_max,
#         em_eff_fwd_range + 0.1,
#         fc_eff_max,
#         # veh_dict['pt_type']['HybridElectricVehicle']['fc'],
#     ])
# )
# assert list(param1_perturb.values()) != list(baseline_errors.values())
# param2_perturb = cal_mod_obj.get_errors(
#     cal_mod_obj.update_params([
#         em_eff_fwd_max,
#         em_eff_fwd_range,
#         fc_eff_max - 0.15,
#         # veh_dict['pt_type']['HybridElectricVehicle']['fc'],
#     ])
# )
# assert list(param2_perturb.values()) != list(baseline_errors.values())
# print("Success!")

if __name__ == "__main__":
    parser = pymoo_api.get_parser()
    args = parser.parse_args()

    n_processes = args.processes 
    n_max_gen = args.n_max_gen 
    # should be at least as big as n_processes
    pop_size = args.pop_size 
    run_minimize = not (args.skip_minimize)

    print("Starting calibration.")
    algorithm = pymoo_api.NSGA2(
        # size of each population
        pop_size=pop_size,
        # LatinHyperCube sampling seems to be more effective than the default
        # random sampling
        sampling=pymoo_api.LHS(),
    )
    termination = pymoo_api.DMOT(
        # max number of generations, default of 10 is very small
        n_max_gen=n_max_gen,
        # evaluate tolerance over this interval of generations every
        period=5,
        # parameter variation tolerance
        xtol=args.xtol,
        # objective variation tolerance
        ftol=args.ftol
    )

    if n_processes == 1:
        print("Running serial evaluation.")
        # series evaluation
        # Setup calibration problem
        cal_prob = pymoo_api.CalibrationProblem(
            mod_obj=cal_mod_obj,
        )
        
        res, res_df = pymoo_api.run_minimize(
            problem=cal_prob,
            algorithm=algorithm,
            termination=termination,
            save_path=save_path,
        )
    else:
        print(f"Running parallel evaluation with n_processes: {n_processes}.")
        assert n_processes > 1
        # parallel evaluation
        import multiprocessing

        with multiprocessing.Pool(n_processes) as pool:
            problem = pymoo_api.CalibrationProblem(
                mod_obj=cal_mod_obj,
                elementwise_runner=StarmapParallelization(pool.starmap),
            )
            res, res_df = pymoo_api.run_minimize(
                problem=problem,
                algorithm=algorithm,
                termination=termination,
                save_path=save_path,
            )

