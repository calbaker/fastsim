"""
Calibration script for 2021_Hyundai_Sonata_Hybrid_Blue    
"""
# # TODO Calibration Tasks
# - [x] put vehicle test data in a sharepoint folder, grant Robin access, and paste link in here
# - [ ] develop means of skewing curves via setter or similar
# - [ ] show what signals should be use for objectives
#     - [x] and how to access them in code
# - [ ] have Robin flesh out and start running calibration
# - [ ] play with things in the meantime

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
sns.set()

# TODO: Kyle or Robin:
# - [ ] in the `./f3-vehicles`, reduce all ~100 element arrays to just the ~10
#       element arrays, and make sure linear interpolation is used
# - [ ] make sure temp- and current/c-rate-dependent battery efficiency interp is being used
veh = fsim.Vehicle.from_file(Path(__file__).parent / "f3-vehicles/2021_Hyundai_Sonata_Hybrid_Blue.yaml")

# Obtain the data from
# https://nrel.sharepoint.com/:f:/r/sites/EEMSCoreModelingandDecisionSupport2022-2024/Shared%20Documents/FASTSim/DynoTestData?csf=1&web=1&e=F4FEBp
# and then copy it to the local folder below
cyc_folder_path = Path(__file__) / "dyno_test_data/2021 Hyundai Sonata Hybrid/Extended Datasets"
assert cyc_folder_path.exists()
cyc_files = [
    # TODO: Someone needs to populate this list after careful review, considering:
    # - HVAC should be on
    # - wide range of initial and ambient temperatures
    # - good signal quality -- somewhat subjective
    # cyc_folder_path / "some_file.txt",
]
assert len(cyc_files) > 0

# TODO: use random selection to retain ~70% of cycles for calibration, and
# reserve the remaining for validation
cyc_files_for_cal = [
    # TOOD: populate this somehow -- e.g. random split of `cyc_files`
]
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
def new_peak_res_eff (sd_dict, new_peak_eff):
    # TODO: fix this dict path because Chad typed it from memory
    sd_dict['veh']['pt_type']['HybridElectricVehicle']['res']['peak_eff'] = new_peak_eff
    # TODO: check that `sd_dict` is mutably modified outside the scope of this function

cal_mod_obj = fsim.pymoo_api.ModelObjectives(
    models = sds_for_cal,
    dfs = dfs_for_cal,
    obj_fns=(
        (
            lambda sd_dict: np.array(sd_dict['veh']['pt_type']['HybridElectricVehicle']['res']['history']['soc']),
            lambda df: df['TODO: find signal for test data soc']
        ),
        # TODO: add objectives for:
        # - engine fuel usage -- 
        # - battery temperature
        # - engine temperature
        # - cabin temperature
    ),
    param_fns=(
        new_peak_res_eff,
        # TODO: make sure this has functions for modifying
        # - battery peak efficiency
        # - engine peak efficiency
        # - motor peak efficiency and efficiency range
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
        # - engine thermal
        #     - thermal mass
        #     - convection to ambient when stopped
        #     - diameter
    ),
    # must match order and length of `params_fns`
    bounds=(
        (0.85, 0.99,),
    ),
    
)

# Setup calibration problem
cal_prob = fsim.pymoo_api.CalibrationProblem()
