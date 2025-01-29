import pandas as pd

from cal_hev import cal_mod_obj, val_mod_obj, save_path

res_df = pd.read_csv(save_path / "pymoo_res_df.csv")
res_df['euclidean'] = (
    res_df.iloc[:, len(cal_mod_obj.param_fns):] ** 2).sum(1).pow(1/2)
best_row = res_df["euclidean"].argmin()
best_df = res_df.iloc[best_row, :]
param_vals = res_df.iloc[best_row, : len(cal_mod_obj.param_fns)].to_numpy()


