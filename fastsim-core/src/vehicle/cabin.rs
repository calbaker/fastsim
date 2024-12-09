use super::*;
// TODO: add parameters and/or cabin model variant for solar heat load

/// Options for handling cabin thermal model
#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, IsVariant)]
pub enum CabinOption {
    /// Basic single thermal capacitance cabin thermal model, including HVAC
    /// system and controls
    LumpedCabin(Box<LumpedCabin>),
    /// Cabin with interior and shell capacitances
    LumpedCabinWithShell,
    /// no cabin thermal model
    #[default]
    None,
}
impl Init for CabinOption {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::LumpedCabin(scc) => scc.init()?,
            Self::LumpedCabinWithShell => {
                todo!()
            }
            Self::None => {}
        }
        Ok(())
    }
}
impl SerdeAPI for CabinOption {}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// Basic single thermal capacitance cabin thermal model, including HVAC
/// system and controls
pub struct LumpedCabin {
    /// Inverse of cabin shell thermal resistance
    pub cab_htc_to_amb: si::HeatTransferCoeff,
    /// parameter for heat transfer coeff from cabin to ambient during
    /// vehicle stop
    pub cab_htc_to_amb_stop: si::HeatTransferCoeff,
    /// cabin thermal capacitance
    pub heat_capacitance: si::HeatCapacity,
    /// cabin length, modeled as a flat plate
    pub length: si::Length,
    /// cabin width, modeled as a flat plate
    pub width: si::Length,
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: LumpedCabinState,
    #[serde(default, skip_serializing_if = "LumpedCabinStateHistoryVec::is_empty")]
    pub history: LumpedCabinStateHistoryVec,
    // TODO: add `save_interval` and associated method
}

impl SerdeAPI for LumpedCabin {}
impl Init for LumpedCabin {}
impl SetCumulative for LumpedCabin {
    fn set_cumulative(&mut self, dt: si::Time) {
        self.state.set_cumulative(dt);
    }
}

impl LumpedCabin {
    /// Solve temperatures, HVAC powers, and cumulative energies of cabin and HVAC system
    /// Arguments:
    /// - `te_amb_air`: ambient air temperature
    /// - `veh_state`: current [VehicleState]
    /// - 'pwr_thermal_from_hvac`: power to cabin from HVAC system
    /// - `dt`: simulation time step size
    pub fn solve(
        &mut self,
        te_amb_air: si::Temperature,
        veh_state: VehicleState,
        pwr_thermal_from_hvac: si::Power,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        self.state.pwr_thermal_from_hvac = pwr_thermal_from_hvac;
        // flat plate model for isothermal, mixed-flow from Incropera and deWitt, Fundamentals of Heat and Mass
        // Transfer, 7th Edition
        let cab_te_film_ext = 0.5 * (self.state.temp + te_amb_air);
        self.state.reynolds_for_plate =
            Air::get_density(Some(cab_te_film_ext), Some(veh_state.elev_curr))
                * veh_state.speed_ach
                * self.length
                / Air::get_dyn_visc(cab_te_film_ext).with_context(|| format_dbg!())?;
        let re_l_crit = 5.0e5 * uc::R; // critical Re for transition to turbulence

        let nu_l_bar: si::Ratio = if self.state.reynolds_for_plate < re_l_crit {
            // equation 7.30
            0.664
                * self.state.reynolds_for_plate.get::<si::ratio>().powf(0.5)
                * Air::get_pr(cab_te_film_ext)
                    .with_context(|| format_dbg!())?
                    .get::<si::ratio>()
                    .powf(1.0 / 3.0)
                * uc::R
        } else {
            // equation 7.38
            let a = 871.0; // equation 7.39
            (0.037 * self.state.reynolds_for_plate.get::<si::ratio>().powf(0.8) - a)
                * Air::get_pr(cab_te_film_ext).with_context(|| format_dbg!())?
        };

        self.state.pwr_thermal_from_amb = if veh_state.speed_ach > 2.0 * uc::MPH {
            let htc_overall_moving: si::HeatTransferCoeff = 1.0
                / (1.0
                    / (nu_l_bar
                        * Air::get_therm_cond(cab_te_film_ext).with_context(|| format_dbg!())?
                        / self.length)
                    + 1.0 / self.cab_htc_to_amb);
            (self.length * self.width) * htc_overall_moving * (te_amb_air - self.state.temp)
        } else {
            (self.length * self.width)
                / (1.0 / self.cab_htc_to_amb_stop + 1.0 / self.cab_htc_to_amb)
                * (te_amb_air - self.state.temp)
        };

        self.state.temp_prev = self.state.temp;
        self.state.temp += (self.state.pwr_thermal_from_hvac + self.state.pwr_thermal_from_amb)
            / self.heat_capacitance
            * dt;
        Ok(())
    }
}

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
pub struct LumpedCabinState {
    /// time step counter
    pub i: u32,
    /// lumped cabin temperature
    pub temp: si::Temperature,
    /// lumped cabin temperature at previous simulation time step
    // TODO: make sure this gets updated
    pub temp_prev: si::Temperature,
    /// Thermal power coming to cabin from HVAC system.  Positive indicates
    /// heating, and negative indicates cooling.
    pub pwr_thermal_from_hvac: si::Power,
    /// Cumulative thermal energy coming to cabin from HVAC system.  Positive indicates
    /// heating, and negative indicates cooling.
    pub energy_thermal_from_hvac: si::Energy,
    /// Thermal power coming to cabin from ambient air.  Positive indicates
    /// heating, and negative indicates cooling.
    pub pwr_thermal_from_amb: si::Power,
    /// Cumulative thermal energy coming to cabin from ambient air.  Positive indicates
    /// heating, and negative indicates cooling.
    pub energy_thermal_from_amb: si::Energy,
    /// Reynolds number for flow over cabin, treating cabin as a flat plate
    pub reynolds_for_plate: si::Ratio,
}

impl Init for LumpedCabinState {}
impl SerdeAPI for LumpedCabinState {}
