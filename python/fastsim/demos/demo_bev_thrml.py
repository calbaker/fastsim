# %%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import seaborn as sns
from pathlib import Path
import time
import json
import os
from typing import Tuple
import fastsim as fsim

sns.set_theme()

from plot_utils  import *

# if environment var `SHOW_PLOTS=false` is set, no plots are shown
SHOW_PLOTS = os.environ.get("SHOW_PLOTS", "true").lower() == "true"     
# if environment var `SAVE_FIGS=true` is set, save plots
SAVE_FIGS = os.environ.get("SAVE_FIGS", "false").lower() == "true"

# `fastsim3` -- load vehicle and cycle, build simulation, and run 
# %%

# load 2020 Chevrolet Bolt BEV from file
veh = fsim.Vehicle.from_file(
    fsim.package_root() / "../../cal_and_val/thermal/f3-vehicles/2020 Chevrolet Bolt EV.yaml"
)

# Set `save_interval` at vehicle level -- cascades to all sub-components with time-varying states
fsim.set_param_from_path(veh, "save_interval" , 1)

# load cycle from file
cyc = fsim.Cycle.from_resource("udds.csv")

# instantiate `SimDrive` simulation object
sd = fsim.SimDrive(veh, cyc)

# simulation start time
t0 = time.perf_counter()
# run simulation
sd.walk()
# simulation end time
t1 = time.perf_counter()
t_fsim3_si1 = t1 - t0
print(f"fastsim-3 `sd.walk()` elapsed time with `save_interval` of 1:\n{t_fsim3_si1:.2e} s")
df = sd.to_dataframe()

# # Visualize results

def plot_res_pwr() -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots(4, 1, sharex=True, figsize=figsize_3_stacked)
    plt.suptitle("Reversible Energy Storage Power")

    ax[0].set_prop_cycle(get_paired_cycler())
    ax[0].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.pwr_out_electrical_watts"] / 1e3,
        label="f3 electrical out",
    )
    ax[0].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.array(sd2.ess_kw_out_ach.tolist()),
        label="f2 electrical out",
    )
    ax[0].set_ylabel("RES Power [kW]")
    ax[0].legend()

    ax[1].set_prop_cycle(get_uni_cycler())
    ax[1].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.pwr_out_electrical_watts"] /
            1e3 - np.array(sd2.ess_kw_out_ach.tolist()),
        label="f3 res kw out",
    )
    ax[1].set_ylabel("RES Power\nDelta (f3-f2) [kW]")
    ax[1].legend()

    ax[2].set_prop_cycle(get_paired_cycler())
    ax[2].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.soc"] - (df["veh.pt_type.BatteryElectricVehicle.res.history.soc"][0] - np.array(sd2.soc.tolist())[0]),
        label="f3 soc",
    )
    ax[2].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.array(sd2.soc.tolist()),
        label="f2 soc",
    )
    ax[2].set_ylabel("SOC")
    ax[2].legend()
    
    ax[-1].set_prop_cycle(get_paired_cycler())
    ax[-1].plot(
        df['cyc.time_seconds'],
        df["veh.history.speed_ach_meters_per_second"],
        label="f3",
    )
    ax[-1].plot(
        np.array(sd2.cyc.time_s.tolist()),
        np.array(sd2.mps_ach.tolist()),
        label="f2",
    )
    ax[-1].legend()
    ax[-1].set_xlabel("Time [s]")
    ax[-1].set_ylabel("Ach Speed [m/s]")

    plt.tight_layout()
    if SAVE_FIGS:
        plt.savefig(Path("./plots/res_pwr.svg"))
    plt.show()

    return fig, ax

def plot_res_energy() -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots(4, 1, sharex=True, figsize=figsize_3_stacked)
    plt.suptitle("Reversible Energy Storage Energy")

    ax[0].set_prop_cycle(get_paired_cycler())
    ax[0].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.energy_out_electrical_joules"] / 1e3,
        label="f3 electrical out",
    )
    ax[0].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.cumsum(np.array(sd2.ess_kw_out_ach.tolist()) * np.diff(sd2.cyc.time_s.tolist(), prepend=0)),
        label="f2 electrical out",
    )
    ax[0].set_ylabel("RES Energy [kW]")
    ax[0].legend()

    ax[1].set_prop_cycle(get_uni_cycler())
    ax[1].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.energy_out_electrical_joules"
            ] / 1e3 - np.cumsum(np.array(sd2.ess_kw_out_ach.tolist()) *
            np.diff(sd2.cyc.time_s.tolist(), prepend=0)),
        label="electrical out",
    )
    ax[1].set_ylim(
       -np.max(np.abs(sd.veh.res.history.energy_out_electrical_joules)) * 1e-3 * 0.1,
        np.max(np.abs(sd.veh.res.history.energy_out_electrical_joules)) * 1e-3 * 0.1
    )
    ax[1].set_ylabel("RES Energy\nDelta (f3-f2) [kJ]\n+/- 10% Range")
    ax[1].legend()

    ax[2].set_prop_cycle(get_paired_cycler())
    ax[2].plot(
        df['cyc.time_seconds'],
        df["veh.pt_type.BatteryElectricVehicle.res.history.soc"] - (df["veh.pt_type.BatteryElectricVehicle.res.history.soc"][0] - np.array(sd2.soc.tolist())[0]),
        label="f3 soc",
    )
    ax[2].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.array(sd2.soc.tolist()),
        label="f2 soc",
    )
    ax[2].set_ylabel("SOC")
    ax[2].legend()
    
    ax[-1].set_prop_cycle(get_paired_cycler())
    ax[-1].plot(
        df['cyc.time_seconds'],
        df["veh.history.speed_ach_meters_per_second"],
        label="f3",
    )
    ax[-1].plot(
        np.array(sd2.cyc.time_s.tolist()),
        np.array(sd2.mps_ach.tolist()),
        label="f2",
    )
    ax[-1].legend()
    ax[-1].set_xlabel("Time [s]")
    ax[-1].set_ylabel("Ach Speed [m/s]")

    plt.tight_layout()
    if SAVE_FIGS:
        plt.savefig(Path("./plots/res_energy.svg"))
    plt.show()

    return fig, ax

def plot_road_loads() -> Tuple[Figure, Axes]: 
    fig, ax = plt.subplots(3, 1, sharex=True, figsize=figsize_3_stacked)
    plt.suptitle("Road Loads")

    ax[0].set_prop_cycle(get_paired_cycler())
    ax[0].plot(
        np.array(sd.cyc.time_seconds)[::veh.save_interval],
        np.array(sd.veh.history.pwr_drag_watts) / 1e3,
        label="f3 drag",
    )
    ax[0].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.array(sd2.drag_kw.tolist()),
        label="f2 drag",
    )
    ax[0].plot(
        np.array(sd.cyc.time_seconds)[::veh.save_interval],
        np.array(sd.veh.history.pwr_rr_watts) / 1e3,
        label="f3 rr",
    )
    ax[0].plot(
        np.array(sd2.cyc.time_s.tolist())[::veh.save_interval],
        np.array(sd2.rr_kw.tolist()),
        label="f2 rr",
    )
    ax[0].set_ylabel("Power [kW]")
    ax[0].legend()

    ax[1].set_prop_cycle(get_uni_cycler())
    ax[1].plot(
        np.array(sd.cyc.time_seconds)[::veh.save_interval],
        np.array(sd.veh.history.pwr_drag_watts) /
        1e3 - np.array(sd2.drag_kw.tolist()),
        label="drag",
        linestyle=BASE_LINE_STYLES[0],
    )
    ax[1].plot(
        np.array(sd.cyc.time_seconds)[::veh.save_interval],
        np.array(sd.veh.history.pwr_rr_watts) /
        1e3 - np.array(sd2.rr_kw.tolist()),
        label="rr",
        linestyle=BASE_LINE_STYLES[1],
    )
    # ax[1].text(
    #     500, -0.125, "Drag error is due to more\naccurate air density model .")
    ax[1].set_ylabel("Power\nDelta (f3-f2) [kW]")
    ax[1].legend()

    ax[-1].set_prop_cycle(get_paired_cycler())
    ax[-1].plot(
        np.array(sd.cyc.time_seconds)[::veh.save_interval],
        np.array(sd.veh.history.speed_ach_meters_per_second),
        label="f3",
    )
    ax[-1].plot(
        np.array(sd2.cyc.time_s.tolist()),
        np.array(sd2.mps_ach.tolist()),
        label="f2",
    )
    ax[-1].legend()
    ax[-1].set_xlabel("Time [s]")
    ax[-1].set_ylabel("Ach. Speed [m/s]")

    plt.tight_layout()
    if SAVE_FIGS:
        plt.savefig(Path("./plots/road_loads.svg"))
    plt.show()

    return fig, ax

if SHOW_PLOTS:
    fig, ax = plot_res_pwr() 
    fig, ax = plot_res_energy()
    fig, ax = plot_road_loads()

# %%
# example for how to use set_default_pwr_interp() method for veh.res
res = veh.res
res.set_default_pwr_interp()
fsim.set_param_from_path(veh, 'res', res)
print(veh.res.to_pydict())

# %%
