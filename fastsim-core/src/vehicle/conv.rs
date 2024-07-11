use super::*;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, HistoryMethods)]
/// Conventional vehicle with only a FuelConverter as a power source
pub struct ConventionalVehicle {
    pub fs: FuelStorage,
    #[has_state]
    pub fc: FuelConverter,
    pub(crate) mass: Option<si::Mass>,
    /// Alternator efficiency used to calculate aux mechanical power demand on engine
    pub alt_eff: si::Ratio,
}

impl SerdeAPI for ConventionalVehicle {}
impl Init for ConventionalVehicle {
    fn init(&mut self) -> anyhow::Result<()> {
        self.fc.init().with_context(|| anyhow!(format_dbg!()))?;
        self.fs.init().with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }
}

impl SaveInterval for ConventionalVehicle {
    fn save_interval(&self) -> anyhow::Result<Option<usize>> {
        bail!("`save_interval` is not implemented in ConventionalVehicle")
    }
    fn set_save_interval(&mut self, save_interval: Option<usize>) -> anyhow::Result<()> {
        self.fc.save_interval = save_interval;
        Ok(())
    }
}

impl Powertrain for Box<ConventionalVehicle> {
    fn set_cur_pwr_prop_out_max(&mut self, pwr_aux: si::Power, dt: si::Time) -> anyhow::Result<()> {
        // TODO: account for transmission efficiency in here
        self.fc
            .set_cur_pwr_tract_out_max(pwr_aux / self.alt_eff, dt)
            .with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }

    fn get_cur_pwr_prop_out_max(&self) -> anyhow::Result<(si::Power, si::Power)> {
        Ok((self.fc.state.pwr_prop_max, 0. * uc::W))
    }

    fn solve(
        &mut self,
        pwr_out_req: si::Power,
        pwr_aux: si::Power,
        _veh_state: &VehicleState,
        _enabled: bool,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        // only positive power can come from powertrain.  Revisit this if engine braking model is needed.
        let pwr_out_req = pwr_out_req.max(uc::W * 0.0);
        let enabled = true; // TODO: replace with a stop/start model
        self.fc
            .solve(pwr_out_req, pwr_aux, enabled, dt)
            .with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }

    fn pwr_regen(&self) -> si::Power {
        uc::W * 0.
    }
}

impl Mass for ConventionalVehicle {
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
                        log::warn!(
                            "Derived mass does not match provided mass, setting `{}` consituent mass fields to `None`",
                            stringify!(ConventionalVehicle));
                        self.expunge_mass_fields();
                    }
                }
                Some(new_mass)
            }
            // Set using `derived_mass()`, failing if it returns `None`
            None => Some(derived_mass.with_context(|| {
                format!(
                    "Not all mass fields in `{}` are set and no mass was provided.",
                    stringify!(ConventionalVehicle)
                )
            })?),
        };
        Ok(())
    }

    fn derived_mass(&self) -> anyhow::Result<Option<si::Mass>> {
        let fc_mass = self.fc.mass().with_context(|| anyhow!(format_dbg!()))?;
        let fs_mass = self.fs.mass().with_context(|| anyhow!(format_dbg!()))?;
        match (fc_mass, fs_mass) {
            (Some(fc_mass), Some(fs_mass)) => Ok(Some(fc_mass + fs_mass)),
            (None, None) => Ok(None),
            _ => bail!(
                "`{}` field masses are not consistently set to `Some` or `None`",
                stringify!(ConventionalVehicle)
            ),
        }
    }

    fn expunge_mass_fields(&mut self) {
        self.fc.mass = None;
        self.fs.mass = None;
        self.mass = None;
    }
}
