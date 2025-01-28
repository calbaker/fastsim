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

import fastsim as fsim

# Initialize seaborn plot configuration
sns.set_style()

veh = fsim.Vehicle.from_file(Path(__file__).parent / "f3-vehicles/2020 Chevrolet Bolt EV.yaml")

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
cyc_files_for_cal = [
    # TOOD: populate this somehow -- e.g. random split of `cyc_files`
]
assert len(cyc_files_for_cal) > 0
dfs_for_cal = {}
for cyc_file in cyc_files_for_cal:
    cyc_file: Path
    # `delimiter="\t"` should work for tab separated variables
    dfs_for_cal[cyc_file.stem] = pd.read_csv(cyc_file, delimiter="\t")
cycs_for_cal = {}
for (cyc_file_stem, df) in dfs_for_cal.items():
    cyc_file_stem: str
    df: pd.DataFrame
    cyc_dict = df.to_dict()
    # TODO: be ready to do some massaging of `cyc_dict`, like making sure that
    # keys match expected, purging invalid keys, and massaging data types

    # TODO: make sure this catches ambient temperature
    cycs_for_cal[cyc_file_stem] = fsim.Cycle.from_pydict(cyc_dict)
sds_for_cal = {}
for (cyc_file_stem, cyc) in cycs_for_cal.items():
    cyc_file_stem: str
    cyc: fsim.Cycle
    # TODO: clone veh and set up initial conditions for:
    # - SOC
    # - cabin temp
    # - battery temp if available, else use cabin temp
    # - engine temp for HEV
    # NOTE: maybe change `save_interval` to 5
    sds_for_cal[cyc_file_stem] = fsim.SimDrive(veh, cyc).to_pydict()

# TODO: flesh this out for validation stuff
# cyc_files_for_val = []

# Setup model objectives
## Parameter Functions
def new_em_eff_max(sd_dict, new_eff_peak):
    """
    Set `new_eff_max` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'])
    em.set_eff_peak(new_eff_peak)
    sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'] = em.to_pydict()
    # TODO: check that `sd_dict` is mutably modified outside the scope of this function, e.g. with a debugger

def new_em_eff_range(sd_dict, new_eff_range):
    """
    Set `new_eff_range` in `ElectricMachine`
    """
    em = fsim.ElectricMachine.from_pydict(sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'])
    em.set_eff_range(new_eff_range)
    sd_dict['veh']['pt_type']['BatteryElectricVehicle']['em'] = em.to_pydict()
    # TODO: check that `sd_dict` is mutably modified outside the scope of this function, e.g. with a debugger

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

# Setup calibration problem
cal_prob = fsim.pymoo_api.CalibrationProblem(
    mod_obj=cal_mod_obj,
)
