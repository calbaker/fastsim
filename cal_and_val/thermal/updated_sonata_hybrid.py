#%%
import fastsim as fsim

veh = fsim.Vehicle.from_resource("2021_Hyundai_Sonata_Hybrid_Blue.yaml")

print(veh.to_pydict())

# %%
#pt_type.HybridElectricVehicle.fc.eff_interp_from_pwr_out.Interp1D.x
# fsim.set_param_from_path(veh, "pt_type.fc.eff_interp_from_pwr_out.x", [0.00, 0.02, 0.04, 0.06, 0.08,  0.10,   0.20,   0.40,   0.60,   0.80,   1.00])
# fsim.set_param_from_path(veh, "pt_type.fc.eff_interp_from_pwr_out.x", [0.83, 0.85,    0.87,   0.89,   0.90,   0.91,   0.93,   0.94,   0.94,   0.93,   0.92])

# %%
veh_pydict = veh.to_pydict()
print(veh_pydict)
# eff_interp_from_pwr_out = fsim.Interpolator.from_pydict(veh_pydict['veh']['pt_type']['BatteryElectricVehicle']['fc']['eff_interp_from_pwr_out'])
# eff_interp_from_pwr_out.set_x([0.00, 0.02, 0.04, 0.06, 0.08,  0.10,   0.20,   0.40,   0.60,   0.80,   1.00])
# eff_interp_from_pwr_out.set_f_x([0.83, 0.85,    0.87,   0.89,   0.90,   0.91,   0.93,   0.94,   0.94,   0.93,   0.92])

# eff_interp_from_pwr_out = veh_pydict['pt_type']['HybridElectricVehicle']['fc']['eff_interp_from_pwr_out']['Interp1D']
# eff_interp_from_pwr_out['x'] = [0.00, 0.02, 0.04, 0.06, 0.08,  0.10,   0.20,   0.40,   0.60,   0.80,   1.00]
# eff_interp_from_pwr_out['f_x'] = [0.83, 0.85,    0.87,   0.89,   0.90,   0.91,   0.93,   0.94,   0.94,   0.93,   0.92]
# veh_pydict['pt_type']['HybridElectricVehicle']['fc']['eff_interp_from_pwr_out'] = eff_interp_from_pwr_out

# veh_pydict['pt_type']['HybridElectricVehicle']['fc']['eff_interp_from_pwr_out']['Interp1D']['x'] = [0.00, 0.02, 0.04, 0.06, 0.08,  0.10,   0.20,   0.40,   0.60,   0.80,   1.00]
# veh_pydict['pt_type']['HybridElectricVehicle']['fc']['eff_interp_from_pwr_out']['Interp1D']['f_x'] = [0.83, 0.85,    0.87,   0.89,   0.90,   0.91,   0.93,   0.94,   0.94,   0.93,   0.92]

veh_pydict['pt_type']['HybridElectricVehicle']['em']['eff_interp_achieved']['Interp1D']['x'] = [0.00, 0.02, 0.04, 0.06, 0.08,  0.10,   0.20,   0.40,   0.60,   0.80,   1.00]
veh_pydict['pt_type']['HybridElectricVehicle']['em']['eff_interp_achieved']['Interp1D']['f_x'] = [0.83, 0.85,    0.87,   0.89,   0.90,   0.91,   0.93,   0.94,   0.94,   0.93,   0.92]

veh = fsim.Vehicle.from_pydict(veh_pydict)
veh.to_yaml()