use super::*;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, HistoryMethods)]
#[non_exhaustive]
/// Battery electric vehicle
pub struct BatteryElectricVehicle {
    #[has_state]
    pub res: ReversibleEnergyStorage,
    #[has_state]
    pub em: ElectricMachine,
    #[has_state]
    pub transmission: Transmission,
    pub(crate) mass: Option<si::Mass>,
}

impl SerdeAPI for BatteryElectricVehicle {}
impl Init for BatteryElectricVehicle {
    fn init(&mut self) -> anyhow::Result<()> {
        self.res.init().with_context(|| anyhow!(format_dbg!()))?;
        self.em.init().with_context(|| anyhow!(format_dbg!()))?;
        self.transmission
            .init()
            .with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }
}

impl Mass for BatteryElectricVehicle {
    fn mass(&self) -> anyhow::Result<Option<si::Mass>> {
        let derived_mass = self
            .derived_mass()
            .with_context(|| anyhow!(format_dbg!()))?;
        match (derived_mass, self.mass) {
            (Some(derived_mass), Some(set_mass)) => {
                ensure!(
                    utils::almost_eq_uom(&set_mass, &derived_mass, None),
                    format!(
                        "{}",
                        format_dbg!(utils::almost_eq_uom(&set_mass, &derived_mass, None)),
                    )
                );
                Ok(Some(set_mass))
            }
            _ => Ok(self.mass.or(derived_mass)),
        }
    }

    fn set_mass(
        &mut self,
        new_mass: Option<si::Mass>,
        side_effect: MassSideEffect,
    ) -> anyhow::Result<()> {
        ensure!(
            side_effect == MassSideEffect::None,
            "At the powertrain level, only `MassSideEffect::None` is allowed"
        );
        let derived_mass = self
            .derived_mass()
            .with_context(|| anyhow!(format_dbg!()))?;
        self.mass = match new_mass {
            // Set using provided `new_mass`, setting constituent mass fields to `None` to match if inconsistent
            Some(new_mass) => {
                if let Some(dm) = derived_mass {
                    if dm != new_mass {
                        self.expunge_mass_fields();
                    }
                }
                Some(new_mass)
            }
            // Set using `derived_mass()`, failing if it returns `None`
            None => Some(derived_mass.with_context(|| {
                format!(
                    "Not all mass fields in `{}` are set and no mass was provided.",
                    stringify!(BatteryElectricVehicle)
                )
            })?),
        };
        Ok(())
    }

    fn derived_mass(&self) -> anyhow::Result<Option<si::Mass>> {
        let res_mass = self.res.mass().with_context(|| anyhow!(format_dbg!()))?;
        let em_mass = self.em.mass().with_context(|| anyhow!(format_dbg!()))?;
        let transmission_mass = self
            .transmission
            .mass()
            .with_context(|| anyhow!(format_dbg!()))?;
        match (res_mass, em_mass, transmission_mass) {
            (Some(res_mass), Some(em_mass), Some(transmission_mass)) => {
                Ok(Some(em_mass + res_mass + transmission_mass))
            }
            (None, None, None) => Ok(None),
            _ => bail!(
                "`{}` field masses are not consistently set to `Some` or `None`",
                stringify!(BatteryElectricVehicle)
            ),
        }
    }

    fn expunge_mass_fields(&mut self) {
        self.res.expunge_mass_fields();
        self.em.expunge_mass_fields();
        self.transmission.expunge_mass_fields();
        self.mass = None;
    }
}

impl SaveInterval for BatteryElectricVehicle {
    fn save_interval(&self) -> anyhow::Result<Option<usize>> {
        bail!("`save_interval` is not implemented in BatteryElectricVehicle")
    }
    fn set_save_interval(&mut self, save_interval: Option<usize>) -> anyhow::Result<()> {
        self.res.save_interval = save_interval;
        self.em.save_interval = save_interval;
        self.transmission.save_interval = save_interval;
        Ok(())
    }
}

impl Powertrain for BatteryElectricVehicle {
    fn solve(
        &mut self,
        pwr_out_req: si::Power,
        _veh_state: VehicleState,
        _enabled: bool,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        let pwr_in_transmission = self
            .transmission
            .get_pwr_in_req(pwr_out_req)
            .with_context(|| anyhow!(format_dbg!()))?;
        let pwr_in_em = self
            .em
            .get_pwr_in_req(pwr_in_transmission, dt)
            .with_context(|| anyhow!(format_dbg!()))?;
        self.res
            .solve(pwr_in_em, dt)
            .with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }

    fn solve_thermal(
        &mut self,
        te_amb: si::Temperature,
        _pwr_thrml_fc_to_cab: si::Power,
        _veh_state: &mut VehicleState,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        self.res
            .solve_thermal(te_amb, dt)
            .with_context(|| format_dbg!())?;
        Ok(())
    }

    fn get_curr_pwr_prop_out_max(&self) -> anyhow::Result<(si::Power, si::Power)> {
        Ok((
            self.em.state.pwr_mech_fwd_out_max,
            self.em.state.pwr_mech_bwd_out_max,
        ))
    }

    fn set_curr_pwr_prop_out_max(
        &mut self,
        pwr_aux: si::Power,
        dt: si::Time,
        _veh_state: VehicleState,
    ) -> anyhow::Result<()> {
        // TODO: account for transmission efficiency in here
        // TODO: change these to something other than zero
        let disch_buffer = si::Energy::ZERO;
        let chrg_buffer = si::Energy::ZERO;
        self.res
            .set_curr_pwr_out_max(dt, disch_buffer, chrg_buffer)
            .with_context(|| anyhow!(format_dbg!()))?;

        self.res
            .set_curr_pwr_prop_max(pwr_aux)
            .with_context(|| anyhow!(format_dbg!()))?;
        self.em
            .set_curr_pwr_prop_out_max(
                self.res.state.pwr_prop_max,
                self.res.state.pwr_regen_max,
                dt,
            )
            .with_context(|| anyhow!(format_dbg!()))?;

        Ok(())
    }

    fn pwr_regen(&self) -> si::Power {
        // When `pwr_mech_prop_out` is negative, regen is happening.  First, clip it at 0, and then negate it.
        // see https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=e8f7af5a6e436dd1163fa3c70931d18d
        // for example
        -self.em.state.pwr_mech_prop_out.min(si::Power::ZERO)
    }
}
