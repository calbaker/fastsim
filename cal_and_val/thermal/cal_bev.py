"""
Calibration script for 2020 Chevrolet Bolt EV    
"""

# critical import
from pathlib import Path

# anticipated cricital imports
import numpy as np  # noqa: F401
import matplotlib.pyplot as plt  # noqa: F401
import seaborn as sns
import pandas as pd  # noqa: F401
import polars as pl  # noqa: F401
from typing import List, Dict
from pymoo.core.problem import StarmapParallelization

import fastsim as fsim
from fastsim import pymoo_api

mps_per_mph = 0.447
celsius_to_kelvin_offset = 273.15

# Initialize seaborn plot configuration
sns.set_style()

veh = fsim.Vehicle.from_file(Path(__file__).parent / "f3-vehicles/2020 Chevrolet Bolt EV.yaml")
veh_dict = veh.to_pydict()

sim_params_dict = fsim.SimParams.default().to_pydict()
sim_params_dict["trace_miss_opts"] = 

# Obtain the data from
# https://nrel.sharepoint.com/:f:/r/sites/EEMSCoreModelingandDecisionSupport2022-2024/Shared%20Documents/FASTSim/DynoTestData?csf=1&web=1&e=F4FEBp
# and then copy it to the local folder below
cyc_folder_path = Path(__file__) / "dyno_test_data/2020 Chevrolet Bolt EV/Extended Datasets"
assert cyc_folder_path.exists()

# See 2020_Chevrolet_Bolt_TestSummary_201005.xlsm for cycle-level data
cyc_files = [
    # TODO: check for seat heater usage in cold cycles and account for that in model!
    # 20F (heater maybe on? Col R in test summary), UDDS + HWY + UDDS + US06
    "62009051 Test Data.txt"
    # 20F (heater maybe on? Col R in test summary), US06 + UDDS + HWY + UDDS
    "62009053 Test Data.txt"

    # room temperature (no HVAC), UDDS + HWY + UDDS + US06
    "62009019 Test Data.txt",
    # room temperature (no HVAC), US06 + UDDS + HWY + UDDS
    "62009021 Test Data.txt",

    # TODO: check for solar load (should be around 1 kW / m^2) and implement or this somewhere (`drive_cycle`???)
    # 95F (HVAC on), UDDS + HWY + UDDS
    "62009040 Test Data.txt"
    # 95F (HVAC on), US06
    "62009041 Test Data.txt"
]
assert len(cyc_files) > 0
cyc_files = [cyc_folder_path / cyc_file for cyc_file in cyc_files]

# TODO: use random selection to retain ~70% of cycles for calibration, and
# reserve the remaining for validation
cyc_files_for_cal: List[str] = [
    # TODO: populate this somehow -- e.g. random split of `cyc_files`
]

cyc_files_for_cal: List[Path] = [cyc_file for cyc_file in cyc_files if cyc_file.name in cyc_files_for_cal]
assert len(cyc_files_for_cal) > 0

def df_to_cyc(df: pd.DataFrame) -> fsim.Cycle:
    # filter out "before" time
    df = df[df["Time[s]_RawFacilities"] >= 0.0]
    assert len(df) > 10
    cyc_dict = {
        "time_seconds": df["Time[s]_RawFacilities"].to_list(),
        "speed_meters_per_second": (df["Dyno_Spd[mph]"] * mps_per_mph).to_list(),
        "temp_amb_air_kelvin": (df["Cell_Temp[C]"] + celsius_to_kelvin_offset).to_list(),
        # TODO: pipe solar load from `Cycle` into cabin thermal model
        # TODO: use something (e.g. regex) to determine solar load
        # see column J comments in 2021_Hyundai_Sonata_Hybrid_TestSummary_2022-03-01_D3.xlsx
        # "pwr_solar_load_watts": df[],
    }
    return fsim.Cycle.from_pydict(cyc_dict, skip_init=False)

def veh_init(cyc_file_stem: str, dfs: Dict[str, pd.DataFrame]) -> fsim.Vehicle:
    # initialize SOC
    # TODO: figure out if `HVBatt_SOC_CAN4__per` is the correct column within the dyno data
    veh_dict['pt_type']['BatteryElectricVehicle']['res']['state']['soc'] = \
        dfs[cyc_file_stem]["HVBatt_SOC_CAN4__per"][0]
    # initialize cabin temp
    veh_dict['cabin']['LumpedCabin']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["Cabin_Temp[C]"][0] + celsius_to_kelvin_offset
    # initialize battery temperature to match cabin temperature because battery
    # temperature is not available in test data
    # Also, battery temperature has no effect in the HEV because efficiency data
    # does not go below 23*C and there is no active thermal management
    veh_dict['pt_type']['BatteryElectricVehicle']['res']['thrml']['RESLumpedThermal']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["Cabin_Temp[C]"][0] + celsius_to_kelvin_offset
    # initialize engine temperature
    veh_dict['pt_type']['BatteryElectricVehicle']['fc']['thrml']['FuelConverterThermal']['state']['temperature_kelvin'] = \
        dfs[cyc_file_stem]["engine_coolant_temp_PCAN__C"][0] + celsius_to_kelvin_offset
    return fsim.Vehicle.from_pydict(veh_dict)

dfs_for_cal: Dict[str, pd.DataFrame] = {
    # `delimiter="\t"` should work for tab separated variables
    cyc_file.stem: pd.read_csv(cyc_file, delimiter="\t") for cyc_file in cyc_files_for_cal
}
for cyc_file in cyc_files_for_cal:
    cyc_file: Path
    # `delimiter="\t"` should work for tab separated variables
    dfs_for_cal[cyc_file.stem] = pd.read_csv(cyc_file, delimiter="\t")
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
    sds_for_cal[cyc_file_stem] = fsim.SimDrive(veh, cyc).to_pydict()

cyc_files_for_val: List[Path] = list(set(cyc_files) - set(cyc_files_for_cal))
assert len(cyc_files_for_val) > 0

dfs_for_val: Dict[str, pd.DataFrame] = {
    # `delimiter="\t"` should work for tab separated variables
    cyc_file.stem: pd.read_csv(cyc_file, delimiter="\t") for cyc_file in cyc_files_for_val
}

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
    sds_for_val[cyc_file_stem] = fsim.SimDrive(veh, cyc).to_pydict()

# Setup model objectives
## Parameter Functions
def new_em_eff_max(sd_dict, new_eff_max):
    """
    Set `new_eff_max` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'])
    em.__eff_fwd_max = new_eff_max
    sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'] = em.to_pydict()
    

def new_em_eff_range(sd_dict, new_eff_range):
    """
    Set `new_eff_range` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'])
    em.__eff_fwd_range = new_eff_range
    sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'] = em.to_pydict()
    

## Model Objectives
cal_mod_obj = fsim.pymoo_api.ModelObjectives(
    models = sds_for_cal,
    dfs = dfs_for_cal,
    obj_fns=(
        (
            lambda sd_dict: np.array(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['res']['history']['soc']),
            lambda df: df['HVBatt_SOC_CAN4__per']
        ),
        # TODO: add objectives for:
        # - battery temperature
        (
            lambda sd_dict: np.array(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['res']['thermal']['RESLumpedThermal']['history']['temperature_kelvin']),
            # HVBatt_cell_temp_1_CAN3__C (or average of temps?) or HVBatt_pack_average_temp_HPCM2__C?
            lambda df: df['HVBatt_pack_average_temp_HPCM2__C']
        ),
        # - cabin temperature
        # - HVAC power, if available
    ),
    param_fns=(
        new_em_eff_max,
        new_em_eff_range,
        # TODO: make sure this has functions for modifying
        # - HVAC PID controls for cabin (not for battery because Sonata has
        #   passive thermal management, but make sure to do battery thermal
        #   controls for BEV)
        # - battery thermal
        #     - thermal mass
        #     - convection to ambient
        #     - convection to cabin
        # - cabin thermal
        #     - thermal mass
        #     - length
        #     - htc to amb when stopped
        #     - set width from vehicle specs -- no need to calibrate
    ),
    # must match order and length of `params_fns`
    bounds=(
        (0.80, 0.99),
        (0.1, 0.6),
    ),
    
)


# verify that model responds to input parameter changes by individually perturbing parameters
baseline_errors = cal_mod_obj.get_errors(
    cal_mod_obj.update_params([
        fsim.ElectricMachine.from_pydict(veh_dict['pt_type']['BatteryElectricVehicle']['em']).eff_fwd_max,
        fsim.ElectricMachine.from_pydict(veh_dict['pt_type']['BatteryElectricVehicle']['em']).eff_fwd_range,
    ])
)
param0_perturb = cal_mod_obj.get_errors(
    cal_mod_obj.update_params([0.90 + 0.5, 0.3])
)
assert list(param0_perturb.values()) != list(baseline_errors.values())
param1_perturb = cal_mod_obj.get_errors(
    cal_mod_obj.update_params([0.90, 0.3 + 0.1])
)
assert list(param1_perturb.values()) != list(baseline_errors.values())

if __name__ == "__main__":
    parser = fsim.cal.get_parser(
        # Defaults are set low to allow for fast run time during testing.  For a good
        # optimization, set this much higher.
        def_save_path=None,
    )
    args = parser.parse_args()

    n_processes = args.processes 
    n_max_gen = args.n_max_gen 
    # should be at least as big as n_processes
    pop_size = args.pop_size 
    run_minimize = not (args.skip_minimize)
    if args.save_path is not None:
        save_path = Path(args.save_path) 
        save_path.mkdir(exist_ok=True)
    else:
        save_path = None

    print("Starting calibration.")
    algorithm = fsim.calibration.NSGA2(
        # size of each population
        pop_size=pop_size,
        # LatinHyperCube sampling seems to be more effective than the default
        # random sampling
        sampling=fsim.calibration.LHS(),
    )
    termination = fsim.calibration.DMOT(
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
            problem = fsim.calibration.CalibrationProblem(
                mod_obj=cal_mod_obj,
                elementwise_runner=StarmapParallelization(pool.starmap),
            )
            res, res_df = pymoo_api.run_minimize(
                problem=problem,
                algorithm=algorithm,
                termination=termination,
                save_path=save_path,
            )


