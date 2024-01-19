use super::*;
// use crate::imports::*;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, HistoryMethods, SerdeAPI)]
/// Battery electric vehicle
pub struct BatteryElectricVehicle {
    #[has_state]
    pub res: ReversibleEnergyStorage,
    #[has_state]
    pub e_machine: ElectricMachine,
}

impl BatteryElectricVehicle {
    pub fn new(
        reversible_energy_storage: ReversibleEnergyStorage,
        e_machine: ElectricMachine,
    ) -> Self {
        BatteryElectricVehicle {
            res: reversible_energy_storage,
            e_machine,
        }
    }
}
impl Powertrain for Box<BatteryElectricVehicle> {
    /// Solve energy consumption for the current power output required
    /// Arguments:
    /// - pwr_out_req: tractive power required
    /// - dt: time step size
    fn solve_powertrain(
        &mut self,
        pwr_out_req: si::Power,
        pwr_aux: si::Power,
        dt: si::Time,
        assert_limits: bool,
    ) -> anyhow::Result<()> {
        self.e_machine.set_pwr_in_req(pwr_out_req, dt)?;
        if self.e_machine.state.pwr_elec_prop_in > si::Power::ZERO {
            // positive traction
            self.res.solve_energy_consumption(
                self.e_machine.state.pwr_elec_prop_in,
                pwr_aux,
                dt,
            )?;
        } else {
            // negative traction
            self.res.solve_energy_consumption(
                self.e_machine.state.pwr_elec_prop_in,
                // limit aux power to whatever is actually available
                // TODO: add more detail/nuance to this
                pwr_aux
                    // whatever power is available from regen plus normal
                    .min(self.res.state.pwr_prop_out_max - self.e_machine.state.pwr_elec_prop_in)
                    .max(si::Power::ZERO),
                dt,
            )?;
        }
        Ok(())
    }

    fn get_pwr_out_max(&mut self, dt: si::Time) -> anyhow::Result<si::Power> {
        todo!();
        // self.res.set_cur_pwr_out_max(pwr_aux.unwrap(), None, None)?;
        // self.e_machine
        //     .set_cur_pwr_max_out(self.res.state.pwr_prop_out_max, None)?;
        // self.e_machine
        //     .set_cur_pwr_regen_max(self.res.state.pwr_regen_out_max)?;

        // // power rate is never limiting in BEL, but assuming dt will be same
        // // in next time step, we can synthesize a rate
        // self.e_machine.set_pwr_rate_out_max(
        //     (self.e_machine.state.pwr_mech_out_max - self.e_machine.state.pwr_mech_prop_out) / dt,
        // );
    }
}
