# Script for validation against fastsim-2

# %%
from pathlib import Path
import fastsim

for veh_file in Path("./f2-vehicles").iterdir():
    print(veh_file)
    fastsim.Vehicle.from_file(veh_file)

# %%
