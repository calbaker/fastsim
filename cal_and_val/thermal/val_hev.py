import pandas as pd
import matplotlib.pyplot as plt

from cal_hev import cal_mod_obj, val_mod_obj, save_path, time_column, mps_per_mph, speed_column

res_df = pd.read_csv(save_path / "pymoo_res_df.csv")
res_df['euclidean'] = (
    res_df.iloc[:, len(cal_mod_obj.param_fns):] ** 2).sum(1).pow(1/2)
best_row = res_df["euclidean"].argmin()
best_df = res_df.iloc[best_row, :]
param_vals = res_df.iloc[best_row, : len(cal_mod_obj.param_fns)].to_numpy()

# getting the solved models
(errors_cal, sds_cal) = cal_mod_obj.get_errors(
    sim_drives=cal_mod_obj.update_params(param_vals),
    return_mods=True,
)
# (errors_val, sds_val) = val_mod_obj.get_errors(
#     sim_drives=val_mod_obj.update_params(param_vals),
#     return_mods=True,
# )

# plotting
plot_save_path = save_path / "plots"
plot_save_path.mkdir(exist_ok=True)

for ((key, df_cal), (sd_key, sd_cal)) in zip(cal_mod_obj.dfs.items(), sds_cal.items()):
    if not isinstance(sd_cal, dict):
        print(f"skipping {key}")
        continue
    assert key == sd_key
    for obj_fn in cal_mod_obj.obj_fns:
        fig, ax = plt.subplots(2, 1, sharex=True)
        fig.suptitle(key)
        ax[0].plot(
            sd_cal['veh']['history']['time_seconds'],
            obj_fn[0](sd_cal),
            label='mod',
        )
        ax[0].plot(
            df_cal[time_column],
            obj_fn[1](df_cal),
            label='exp',
        )
        ax[0].legend()
        ax[0].set_ylabel(obj_fn[0].__name__)

        ax[1].plot(
            sd_cal['veh']['history']['time_seconds'],
            sd_cal['veh']['history']['speed_ach_meters_per_second'],
            label='mod',
        )
        ax[1].plot(
            df_cal[time_column],
            df_cal[speed_column] * mps_per_mph,
            label='exp',
        )
        ax[1].legend()
        ax[1].set_ylabel("Speed [m/s]")
        plt.savefig(plot_save_path / f"{key}.svg")
