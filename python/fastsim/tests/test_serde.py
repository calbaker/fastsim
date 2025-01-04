import numpy as np
from pathlib import Path
import time
import json
import msgpack
import os
from typing import Tuple
import fastsim as fsim

# load 2012 Ford Fusion from file
veh = fsim.Vehicle.from_resource("2022_Renault_Zoe_ZE50_R135.yaml")

# Set `save_interval` at vehicle level -- cascades to all sub-components with time-varying states
fsim.set_param_from_path(veh, "save_interval" , 1)

# load cycle from file
cyc = fsim.Cycle.from_resource("udds.csv")

# instantiate SimDrive 
sd = fsim.SimDrive(veh, cyc)

# simulation start time
# run simulation
sd.walk()
# simulation end time

# TODO: make performmance benchmarks comparing full circle (de)serialization
# with different formats (e.g. message pack, json, yaml, ...)

def test_msg_pack():
    sd_dict = msgpack.loads(sd.to_msg_pack())
    assert sd.to_msg_pack() == msgpack.packb(sd_dict)
    fsim.SimDrive.from_msg_pack(msgpack.packb(sd_dict))
