use fmt::Display;

use crate::prelude::{ElectricMachineState, FuelConverterState};

use super::{vehicle_model::VehicleState, *};

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize, HistoryMethods)]
/// Hybrid vehicle with both engine and reversible energy storage (aka battery)
/// This type of vehicle is not likely to be widely prevalent due to modularity of consists.
pub struct HybridElectricVehicle {
    #[has_state]
    pub res: ReversibleEnergyStorage,
    pub fs: FuelStorage,
    #[has_state]
    pub fc: FuelConverter,
    #[has_state]
    pub em: ElectricMachine,
    /// control strategy for distributing power demand between `fc` and `res`
    #[serde(default)]
    pub pt_cntrl: HEVPowertrainControls,
    /// control strategy for distributing aux power demand between `fc` and `res`
    #[serde(default)]
    pub aux_cntrl: HEVAuxControls,
    /// hybrid powertrain mass
    pub(crate) mass: Option<si::Mass>,
    #[serde(default)]
    pub sim_params: HEVSimulationParams,
    /// Number of `walk` iterations required to achieve SOC balance (i.e. SOC
    /// ends at same starting value, ensuring no net [ReversibleEnergyStorage] usage)
    #[serde(default)]
    pub soc_bal_iters: u32,
    /// field for tracking current state
    #[serde(default)]
    #[serde(skip_serializing_if = "EqDefault::eq_default")]
    pub state: HEVState,
    /// vector of [Self::state]
    #[serde(default)]
    #[serde(skip_serializing_if = "HEVStateHistoryVec::is_empty")]
    pub history: HEVStateHistoryVec,
}

impl SaveInterval for HybridElectricVehicle {
    fn save_interval(&self) -> anyhow::Result<Option<usize>> {
        bail!("`save_interval` is not implemented in HybridElectricVehicle")
    }
    fn set_save_interval(&mut self, save_interval: Option<usize>) -> anyhow::Result<()> {
        self.res.save_interval = save_interval;
        self.em.save_interval = save_interval;
        Ok(())
    }
}

impl Init for HybridElectricVehicle {
    fn init(&mut self) -> anyhow::Result<()> {
        self.fc.init().with_context(|| anyhow!(format_dbg!()))?;
        self.res.init().with_context(|| anyhow!(format_dbg!()))?;
        self.em.init().with_context(|| anyhow!(format_dbg!()))?;
        self.pt_cntrl
            .init()
            .with_context(|| anyhow!(format_dbg!()))?;
        self.state.init().with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }
}

impl Powertrain for Box<HybridElectricVehicle> {
    fn set_curr_pwr_prop_out_max(
        &mut self,
        pwr_aux: si::Power,
        dt: si::Time,
        veh_state: &VehicleState,
    ) -> anyhow::Result<()> {
        // TODO: account for transmission efficiency in here
        self.state.fc_on_causes.clear();
        match &self.pt_cntrl {
            HEVPowertrainControls::Fastsim2(rgwb) => {
                if self.fc.state.fc_on && self.fc.state.time_on
                    < rgwb.fc_min_time_on.with_context(|| {
                    anyhow!(
                        "{}\n Expected `ResGreedyWithBuffers::init` to have been called beforehand.",
                        format_dbg!()
                    )
                })? {
                    self.state.fc_on_causes.push(FCOnCause::OnTimeTooShort)
                }
            },
            HEVPowertrainControls::RESGreedyWithDynamicBuffers => {
                todo!()
            }
        };
        self.fc
            .set_curr_pwr_out_max(dt)
            .with_context(|| anyhow!(format_dbg!()))?;
        let disch_buffer = match &self.pt_cntrl {
            HEVPowertrainControls::Fastsim2(rgwb) => {
                if veh_state.speed_ach
                    < rgwb.speed_soc_buffer_for_accel_cutoff.with_context(|| {
                        anyhow!(
                        "{}\nExpected `ResGreedyWithBuffers::init` to have been called beforehand.",
                        format_dbg!()
                    )
                    })?
                {
                    rgwb.soc_frac_buffer_for_accel.with_context(|| {
                        anyhow!(
                        "{}\nExpected `ResGreedyWithBuffers::init` to have been called beforehand.",
                        format_dbg!()
                    )
                    })? * self.res.energy_capacity_usable()
                } else {
                    si::Energy::ZERO
                }
            }
            HEVPowertrainControls::RESGreedyWithDynamicBuffers => {
                todo!()
            }
        };
        // TODO: change this to something other than zero
        let chrg_buffer = si::Energy::ZERO;
        self.res
            .set_curr_pwr_out_max(dt, disch_buffer, chrg_buffer)
            .with_context(|| anyhow!(format_dbg!()))?;
        let (pwr_aux_res, pwr_aux_fc) = {
            match self.aux_cntrl {
                HEVAuxControls::AuxOnResPriority => {
                    if pwr_aux <= self.res.state.pwr_disch_max {
                        (pwr_aux, si::Power::ZERO)
                    } else {
                        (si::Power::ZERO, pwr_aux)
                    }
                }
                HEVAuxControls::AuxOnFcPriority => (si::Power::ZERO, pwr_aux),
            }
        };
        if pwr_aux_fc > si::Power::ZERO {
            self.state.fc_on_causes.push(FCOnCause::AuxPowerDemand);
        }
        self.fc
            .set_curr_pwr_prop_max(pwr_aux_fc)
            .with_context(|| anyhow!(format_dbg!()))?;
        self.res
            .set_curr_pwr_prop_max(pwr_aux_res)
            .with_context(|| anyhow!(format_dbg!()))?;
        self.em
            .set_curr_pwr_prop_out_max(
                // TODO: add means of controlling whether fc can provide power to em and also how much
                // Try out a 'power out type' enum field on the fuel converter with variants for mechanical and electrical
                self.res.state.pwr_prop_max,
                self.res.state.pwr_regen_max,
                dt,
            )
            .with_context(|| anyhow!(format_dbg!()))?;
        Ok(())
    }

    fn get_curr_pwr_prop_out_max(&self) -> anyhow::Result<(si::Power, si::Power)> {
        Ok((
            self.em.state.pwr_mech_fwd_out_max + self.fc.state.pwr_prop_max,
            self.em.state.pwr_mech_bwd_out_max,
        ))
    }

    fn solve(
        &mut self,
        pwr_out_req: si::Power,
        _veh_state: &VehicleState,
        _enabled: bool,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        let (fc_pwr_out_req, em_pwr_out_req) = self.pt_cntrl.get_pwr_fc_and_em(
            pwr_out_req,
            &mut self.state,
            &mut self.fc.state,
            &self.em.state,
        )?;
        // TODO: add a stop/start model
        // TODO: figure out fancier way to handle apportionment of `pwr_aux` between `fc` and `res`

        let fc_on: bool = !self.state.fc_on_causes.0.is_empty(); // placeholder, delete on merge conflict

        self.fc
            .solve(fc_pwr_out_req, fc_on, dt)
            .with_context(|| format_dbg!())?;
        let res_pwr_out_req = self
            .em
            .get_pwr_in_req(em_pwr_out_req, dt)
            .with_context(|| format_dbg!())?;
        self.res
            .solve(res_pwr_out_req, dt)
            .with_context(|| format_dbg!())?;
        Ok(())
    }

    fn pwr_regen(&self) -> si::Power {
        // When `pwr_mech_prop_out` is negative, regen is happening.  First, clip it at 0, and then negate it.
        // see https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=e8f7af5a6e436dd1163fa3c70931d18d
        // for example
        -self.em.state.pwr_mech_prop_out.min(0. * uc::W)
    }
}

impl Mass for HybridElectricVehicle {
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
                    stringify!(HybridElectricVehicle)
                )
            })?),
        };
        Ok(())
    }

    fn derived_mass(&self) -> anyhow::Result<Option<si::Mass>> {
        let fc_mass = self.fc.mass().with_context(|| anyhow!(format_dbg!()))?;
        let fs_mass = self.fs.mass().with_context(|| anyhow!(format_dbg!()))?;
        let res_mass = self.res.mass().with_context(|| anyhow!(format_dbg!()))?;
        let em_mass = self.em.mass().with_context(|| anyhow!(format_dbg!()))?;
        match (fc_mass, fs_mass, res_mass, em_mass) {
            (Some(fc_mass), Some(fs_mass), Some(res_mass), Some(em_mass)) => {
                Ok(Some(fc_mass + fs_mass + em_mass + res_mass))
            }
            (None, None, None, None) => Ok(None),
            _ => bail!(
                "`{}` field masses are not consistently set to `Some` or `None`",
                stringify!(HybridElectricVehicle)
            ),
        }
    }

    fn expunge_mass_fields(&mut self) {
        self.fc.expunge_mass_fields();
        self.fs.expunge_mass_fields();
        self.res.expunge_mass_fields();
        self.em.expunge_mass_fields();
        self.mass = None;
    }
}

#[fastsim_api]
#[derive(Clone, Debug, Default, Deserialize, PartialEq)]
pub struct FCOnCauses(Vec<FCOnCause>);
impl Init for FCOnCauses {}
impl SerdeAPI for FCOnCauses {}
impl FCOnCauses {
    fn clear(&mut self) {
        self.0.clear();
    }

    #[allow(dead_code)]
    fn pop(&mut self) -> Option<FCOnCause> {
        self.0.pop()
    }

    fn push(&mut self, new: FCOnCause) {
        self.0.push(new)
    }
}

// TODO: figure out why this is not turning in the dataframe but is in teh pydict
#[fastsim_api]
#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative)]
pub struct HEVState {
    /// time step index
    pub i: usize,
    /// Vector of posssible reasons the fc is forced on
    #[api(skip_get, skip_set)]
    pub fc_on_causes: FCOnCauses,
    // TODO: store the soc buffer here or somewhere
}
impl Init for HEVState {}
impl SerdeAPI for HEVState {}

// TODO: implement `Deserialize`

// Custom serialization
impl Serialize for FCOnCauses {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let joined = self
            .0
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<String>>()
            .join(", ");
        serializer.serialize_str(&format!("\"[{}]\"", joined))
    }
}
impl Display for FCOnCauses {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
        // or, alternatively:
        // fmt::Debug::fmt(self, f)
    }
}

#[fastsim_enum_api]
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub enum FCOnCause {
    /// Engine must be on to self heat if thermal model is enabled
    FCTemperatureTooLow,
    /// Engine must be on for high vehicle speed to ensure powertrain can meet
    /// any spikes in power demand
    VehicleSpeedTooHigh,
    /// Engine has not been on long enough (usually 30 s)
    OnTimeTooShort,
    /// Powertrain power demand exceeds motor and/or battery capabilities
    PropulsionPowerDemand,
    /// Aux power demand exceeds battery capability
    AuxPowerDemand,
}
impl SerdeAPI for FCOnCause {}
impl Init for FCOnCause {}
impl fmt::Display for FCOnCause {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
        // or, alternatively:
        // fmt::Debug::fmt(self, f)
    }
}

/// Options for controlling simulation behavior
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub struct HEVSimulationParams {
    /// [ReversibleEnergyStorage] per [FuelConverter]
    pub res_per_fuel_lim: si::Ratio,
    /// Threshold of SOC balancing iterations for triggering warning
    pub soc_balance_iter_warn: u32,
    /// Threshold of SOC balancing iteration for triggering error
    pub soc_balance_iter_err: u32,
    /// Whether to allow iteration to achieve SOC balance
    pub balance_soc: bool,
}

impl Default for HEVSimulationParams {
    fn default() -> Self {
        Self {
            res_per_fuel_lim: uc::R * 0.005,
            soc_balance_iter_warn: 3,
            soc_balance_iter_err: 5,
            balance_soc: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize, Default)]
pub enum HEVAuxControls {
    /// If feasible, use [ReversibleEnergyStorage] to handle aux power demand
    #[default]
    AuxOnResPriority,
    /// If feasible, use [FuelConverter] to handle aux power demand
    AuxOnFcPriority,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub enum HEVPowertrainControls {
    /// Controls that attempt to match fastsim-2
    Fastsim2(RESGreedyWithBuffers),
    /// Controls that have a dynamically updated discharge buffer but are otherwise similar to [Self::Fastsim2]
    RESGreedyWithDynamicBuffers,
}

impl Default for HEVPowertrainControls {
    fn default() -> Self {
        Self::Fastsim2(Default::default())
    }
}

impl Init for HEVPowertrainControls {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::Fastsim2(rgwb) => rgwb.init()?,
            Self::RESGreedyWithDynamicBuffers => {
                todo!()
            }
        }
        Ok(())
    }
}

impl HEVPowertrainControls {
    fn get_pwr_fc_and_em(
        &self,
        pwr_out_req: si::Power,
        hev_state: &mut HEVState,
        fc_state: &mut FuelConverterState,
        em_state: &ElectricMachineState,
    ) -> anyhow::Result<(si::Power, si::Power)> {
        if pwr_out_req >= si::Power::ZERO {
            // positive net power out of the powertrain
            // cannot exceed ElectricMachine max output power. Excess demand will be handled by `fc`
            let em_pwr = pwr_out_req.min(em_state.pwr_mech_fwd_out_max);

            let fc_pwr_frac_demand_forced_on = match &self {
                HEVPowertrainControls::Fastsim2(rgwb) => rgwb
                    .fc_pwr_frac_demand_forced_on
                    .with_context(|| format_dbg!())?,
                HEVPowertrainControls::RESGreedyWithDynamicBuffers => uc::R * f64::NAN,
            };
            let fc_pwr = pwr_out_req - em_pwr;
            // If the motor cannot produce more than the required power times a
            // 0..1 fraction, then the engine should be on
            if em_state.pwr_mech_fwd_out_max < pwr_out_req * fc_pwr_frac_demand_forced_on
                || fc_pwr > si::Power::ZERO
            {
                hev_state
                    .fc_on_causes
                    .push(FCOnCause::PropulsionPowerDemand);
            }

            ensure!(
                fc_pwr >= si::Power::ZERO,
                format_dbg!(fc_pwr >= si::Power::ZERO)
            );
            ensure!(
                pwr_out_req <= em_state.pwr_mech_fwd_out_max + fc_state.pwr_prop_max,
                "{}\n`pwr_out_req`: {} kW\n`em_state.pwr_mech_fwd_out_max`: {} kW",
                format_dbg!(pwr_out_req <= em_state.pwr_mech_fwd_out_max),
                pwr_out_req.get::<si::kilowatt>(),
                em_state.pwr_mech_fwd_out_max.get::<si::kilowatt>()
            );
            Ok((fc_pwr, em_pwr))
        } else {
            // negative net power out of the powertrain -- i.e. positive net
            // power _into_ powertrain, aka regen
            // if `em_pwr` is less than magnitude of `pwr_out_req`, friction brakes can handle excess
            let em_pwr = -em_state.pwr_mech_bwd_out_max.min(-pwr_out_req);
            Ok((0. * uc::W, em_pwr))
        }
    }
}

// TODO: make sure this code has an equivalent in fastsim-3
// if self.veh.no_elec_sys {
//     self.regen_buff_soc[i] = 0.0;
// } else if self.veh.charging_on {
//     self.regen_buff_soc[i] = max(
//         self.veh.max_soc - (self.veh.max_regen_kwh / self.veh.ess_max_kwh),
//         (self.veh.max_soc + self.veh.min_soc) / 2.0,
//     );
// } else {
//     self.regen_buff_soc[i] = max(
//         (self.veh.ess_max_kwh * self.veh.max_soc
//             - 0.5
//                 * self.veh.veh_kg
//                 * (self.cyc.mps[i].powi(2))
//                 * (1.0 / 1_000.0)
//                 * (1.0 / 3_600.0)
//                 * self.veh.mc_peak_eff()
//                 * self.veh.max_regen)
//             / self.veh.ess_max_kwh,
//         self.veh.min_soc,
//     );

/// Container for static controls parameters
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize, Default)]
pub struct RESGreedyWithBuffers {
    /// Fraction of usable battery energy above static minimum SOC that is reserved for acceleration.  Defaults to 0.2.
    pub soc_frac_buffer_for_accel: Option<si::Ratio>,
    /// Speed at which accel buffer becomes inactive.  Defaults to 60 mph.
    // NOTE: for future control strategy, make this a ramp rather than a cutoff
    pub speed_soc_buffer_for_accel_cutoff: Option<si::Velocity>,
    /// Minimum time engine must remain on if it was on during the previous
    /// simulation time step.  defaults to 30 s.
    pub fc_min_time_on: Option<si::Time>,
    /// Speed at which [fuelconverter] is forced on. Defaults to 75 mph.
    //  TODO: make sure this is plumbed up
    pub fc_speed_forced_on: Option<si::Velocity>,
    /// Fraction of total aux and powertrain power demand at which [FuelConverter] is forced on.  Defaults to 0.75.
    // TODO: make sure this is plumbed up
    pub fc_pwr_frac_demand_forced_on: Option<si::Ratio>,
    /// Fraction of available charging capacity to use toward running the engine efficiently. Defaults to 0.
    // TODO: make sure this is plumbed up
    pub frac_res_chrg_for_fc: si::Ratio,
    /// Fraction of available discharging capacity to use toward running the engine efficiently. Defaults to 0.
    // TODO: make sure this is plumbed up
    pub frac_res_dschrg_for_fc: si::Ratio,
    /// Fraction of usable battery energy to reserve for charging when decelerating from high speed
    // TODO: make sure this is plumbed up
    pub soc_frac_buffer_for_regen: Option<si::Ratio>,
}

impl Init for RESGreedyWithBuffers {
    fn init(&mut self) -> anyhow::Result<()> {
        self.soc_frac_buffer_for_accel = self.soc_frac_buffer_for_accel.or(Some(0.2 * uc::R));
        self.speed_soc_buffer_for_accel_cutoff = self
            .speed_soc_buffer_for_accel_cutoff
            .or(Some(60. * uc::MPH));
        self.fc_min_time_on = self.fc_min_time_on.or(Some(uc::S * 30.));
        self.fc_speed_forced_on = self.fc_speed_forced_on.or(Some(uc::MPH * 75.));
        self.fc_pwr_frac_demand_forced_on =
            self.fc_pwr_frac_demand_forced_on.or(Some(uc::R * 0.75));
        Ok(())
    }
}
