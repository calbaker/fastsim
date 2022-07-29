# %%
from pymoo.util.termination.default import MultiObjectiveDefaultTermination as MODT
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np

import fastsim as fsim
import fastsimrust as fsr

# load test data which can be obtained at
# https://www.anl.gov/taps/d3-2012-ford-fusion-v6
possible_trip_dirs = (
    Path().home() / "Documents/DynoTestData/FordFusionTestData/",
)

for trip_dir in possible_trip_dirs:
    if trip_dir.exists():
        break

rho_fuel_kg_per_ml = 0.743e-3
lhv_fuel_btu_per_lbm = 18_344
lbm_per_kg = 2.2
btu_per_kj = 0.948
lhv_fuel_kj_per_kg = lhv_fuel_btu_per_lbm * lbm_per_kg / btu_per_kj

# full data
dfs_raw = dict()
# resampled to 1 Hz
dfs = dict()
for sub in trip_dir.iterdir():
    if sub.is_dir():
        for file in sub.iterdir():
            if file.suffix == ".csv" and "_cs" in file.stem:
                print(f"loading: ", file.resolve())
                dfs_raw[file.stem] = pd.read_csv(file)
                # clip time at zero seconds
                dfs_raw[file.stem] = dfs_raw[file.stem][dfs_raw[file.stem]
                                                        ['Time[s]'] >= 0.0]
                dfs[file.stem] = fsim.resample(
                    dfs_raw[file.stem],
                    rate_vars=('Eng_FuelFlow_Direct[cc/s]')
                )
                dfs[file.stem]['Fuel_Power_Calc[kW]'] = dfs[
                    file.stem]["Eng_FuelFlow_Direct[cc/s]"] * rho_fuel_kg_per_ml * lhv_fuel_kj_per_kg

# %%
# plot the data
show_plots = False

if show_plots:
    for key, df in dfs_raw.items():
        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 6))
        ax[0].set_title(key)
        ax[0].plot(df['Time[s]'], df["CylinderHeadTempC"], label="cyl. head")
        ax[0].plot(df['Time[s]'], df["Cell_Temp[C]"], label="ambient")
        ax[0].set_ylabel("temp [°C]")
        ax[0].legend()
        ax[1].plot(df['Time[s]'], df["Fuel_Power_Calc[kW]"], label="ambient")
        ax[1].set_ylabel("Fuel Power [kW]")
        ax[1].legend()
        ax[-1].plot(df['Time[s]'], df["Dyno_Spd[mph]"])
        ax[-1].set_ylabel("speed [mph]")
        ax[-1].set_xlabel('time [s]')

# %%
# Separate calibration and validation cycles

cal_cyc_patterns = ("49", "56", "73", "60", "69", "77")
dfs_cal = dict()
for key in dfs.keys():
    for pattern in cal_cyc_patterns:
        if pattern in key:
            dfs_cal[key] = dfs[key]

dfs_val_keys = set(dfs.keys()) - set(dfs_cal.keys())
dfs_val = {key: dfs[key] for key in dfs_val_keys}


# %%
# create cycles and sim_drives

veh = fsim.vehicle.Vehicle.from_file("2012_Ford_Fusion.csv").to_rust()

cycs = dict()
cal_sim_drives = dict()
val_sim_drives = dict()
for key in dfs.keys():
    cycs[key] = fsim.cycle.Cycle.from_dict(
        {
            "time_s": dfs[key]["Time[s]"],
            "mps": dfs[key]["Dyno_Spd[mph]"] / fsim.params.MPH_PER_MPS
        }
    ).to_rust()
    sdh = fsr.SimDriveHot(
        cycs[key],
        veh,
        amb_te_deg_c=dfs[key]['Cell_Temp[C]'][0],
        fc_te_deg_c_init=dfs[key]['CylinderHeadTempC'][0])
    sim_params = sdh.sd.sim_params
    sim_params.reset_orphaned()
    # make tolerances big since runs may print lots of warnings before final design is selected
    sim_params.trace_miss_speed_mps_tol = 1e9
    sim_params.trace_miss_time_tol = 1e9
    sim_params.trace_miss_dist_tol = 1e9
    sim_params.energy_audit_error_tol = 1e9
    sim_params.verbose = False
    sd = sdh.sd
    sd.reset_orphaned()
    sd.sim_params = sim_params
    sdh.sd = sd

    if key in list(dfs_cal.keys()):
        cal_sim_drives[key] = sdh
    else:
        assert key in list(dfs_val.keys())
        val_sim_drives[key] = sdh


# %%
objectives = fsim.calibration.ModelErrors(
    sim_drives=cal_sim_drives,
    dfs=dfs_cal,
    objectives=[
        ("sd.fs_kw_out_ach", "Fuel_Power_Calc[kW]"),
        ("history.fc_te_deg_c", "CylinderHeadTempC"),
    ],
    params=[
        "vehthrm.fc_c_kj__k",
        "vehthrm.fc_l",
        "vehthrm.fc_htc_to_amb_stop",
        "vehthrm.fc_coeff_from_comb",
        "vehthrm.fc_exp_offset",
        "vehthrm.fc_exp_lag",
        "vehthrm.fc_exp_minimum",
        "vehthrm.rad_eps",
    ],
    verbose=False
)

# Kyle, when you get this running, remind me to show you how to parallelize it.

problem = fsim.calibration.CalibrationProblem(
    err=objectives,
    param_bounds=[
        (50, 200),
        (0.25, 2),
        (5, 50),
        (1e-5, 1e-3),
        (-10, 30),
        (15, 75),
        (0.25, 0.45),
        (5, 50),
    ],
    # n_max_gen=100,
    # pop_size=12,
)

# %%

b_run_min = False

res_save_path = "pymoo_res_df.csv"

if b_run_min:
    res, res_df = fsim.calibration.run_minimize(
        problem,
        termination=MODT(
            # max number of generations, default of 10 is very small
            n_max_gen=500,
            # evaluate tolerance over this interval of generations every `nth_gen`
            n_last=10,
        ),
    )
    res_df
else:
    res_df = pd.read_csv(res_save_path)

res_df['euclidean'] = (
    res_df.iloc[:, len(objectives.params):] ** 2).sum(1).pow(1/2)
best_row = res_df['euclidean'].argmin()
best_df = res_df.iloc[best_row, :]
param_vals = res_df.iloc[0, :len(objectives.params)].to_numpy()

# %%

# params_and_vals = {
#     'vehthrm.fc_c_kj__k': 125.0,
#     'vehthrm.fc_l': 1.3,
#     'vehthrm.fc_htc_to_amb_stop': 100.0,
#     'vehthrm.fc_coeff_from_comb': 0.00030721481819805005,
#     'vehthrm.fc_exp_offset': -9.438669088889137,
#     'vehthrm.fc_exp_lag': 30.0,
#     'vehthrm.fc_exp_minimum': 0.2500008623533276,
#     'vehthrm.rad_eps': 20
#  }

plot_save_dir = Path("plots")
plot_save_dir.mkdir(exist_ok=True)

# problem.err.update_params(params_and_vals.values())
problem.err.update_params(param_vals)
problem.err.get_errors(plot=True, plot_save_dir=plot_save_dir, plot_perc_err=False)


# %%

# Demonstrate with model showing fuel usage impact

te_amb_deg_c_arr = np.arange(-10, 51)
mpg_arr = np.zeros(len(te_amb_deg_c_arr))

for i, te_amb_deg_c in enumerate(te_amb_deg_c_arr):
    sdh = fsr.SimDriveHot(
        fsim.cycle.Cycle.from_file("udds").to_rust(),
        veh,
        amb_te_deg_c=te_amb_deg_c,
        fc_te_deg_c_init=te_amb_deg_c
    )
    sdh.sim_drive()
    mpg_arr[i] = sdh.sd.mpgge

# %%
fig, ax = plt.subplots()
ax.scatter(te_amb_deg_c_arr, mpg_arr)
ax.set_xlabel("Ambient/Cold Start Temperature [°C]")
ax.set_ylabel("Fuel Economy [mpg]")
ax.set_title("2012 Ford Fusion V6")
plt.tight_layout()
plt.savefig("plots/fe v amb temp.svg")
plt.savefig("plots/fe v amb temp.png")

