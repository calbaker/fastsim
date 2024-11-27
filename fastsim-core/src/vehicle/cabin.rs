use super::*;

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
            Self::LumpedCabinWithShell => {}
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
    /// cabin shell thermal resistance \[m **2 * K / W\]
    pub cab_r_to_amb: f64,
    /// parameter for heat transfer coeff \[W / (m ** 2 * K)\] from cabin to ambient during
    /// vehicle stop
    pub cab_htc_to_amb_stop: f64,
    /// cabin thermal capacitance
    pub heat_capacitance: si::HeatCapacity,
    /// cabin length, modeled as a flat plate
    pub length: si::Length,
    /// cabin width, modeled as a flat plate
    pub width: si::Length,
    /// HVAC model
    pub hvac: HVACSystemForSCC,
    pub state: LumpedCabinState,
    #[serde(default)]
    #[serde(skip_serializing_if = "LumpedCabinStateHistoryVec::is_empty")]
    pub history: LumpedCabinStateHistoryVec,
}

impl SerdeAPI for LumpedCabin {}
impl Init for LumpedCabin {}
impl SetCumulative for LumpedCabin {
    fn set_cumulative(&mut self, dt: si::Time) {
        self.state.set_cumulative(dt);
    }
}

impl LumpedCabin {
    /// Solve temperatures, powers, and cumulative energies of cabin and HVAC system
    /// Arguments:
    /// - `te_amb_air`: ambient air temperature
    /// - `te_fc`: [FuelConverter] temperature, as appropriate for [PowertrainType]
    /// - `dt`: simulation time step size
    fn solve(
        &mut self,
        te_amb_air: si::Temperature,
        te_fc: Option<si::Temperature>,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        if self.state.temp <= self.hvac.te_set + self.hvac.te_deadband
            && self.state.temp >= self.hvac.te_set - self.hvac.te_deadband
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_thermal_from_hvac = si::Power::ZERO;
            self.hvac.state.pwr_i = si::Power::ZERO; // reset to 0.0
            self.hvac.state.pwr_p = si::Power::ZERO;
            self.hvac.state.pwr_d = si::Power::ZERO;
        } else {
            let te_delta_vs_set = self.state.temp - self.hvac.te_set;
            let te_delta_vs_amb: si::Temperature = self.state.temp - te_amb_air;

            self.hvac.state.pwr_p = self.hvac.p * te_delta_vs_set;
            self.hvac.state.pwr_i +=
                (self.hvac.i * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt)
                    .max(self.hvac.pwr_i_max);
            self.hvac.state.pwr_d =
                self.hvac.d * uc::J / uc::KELVIN * ((self.state.temp - self.state.temp_prev) / dt);

            // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
            // cop_ideal is t_h / (t_h - t_c) for heating
            // cop_ideal is t_c / (t_h - t_c) for cooling

            // divide-by-zero protection and realistic limit on COP
            let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                // cabin is cooler than ambient + threshold
                // TODO: make this `5.0` not hardcoded
                self.state.temp / (5.0 * uc::KELVIN)
            } else {
                self.state.temp / te_delta_vs_amb.abs()
            };
            self.hvac.state.cop = cop_ideal * self.hvac.frac_of_ideal_cop;
            assert!(self.hvac.state.cop > 0.0 * uc::R);

            if self.state.temp > self.hvac.te_set + self.hvac.te_deadband {
                // COOLING MODE; cabin is hotter than set point

                if self.hvac.state.pwr_i < si::Power::ZERO {
                    // reset to switch from heating to cooling
                    self.hvac.state.pwr_i = si::Power::ZERO;
                }
                self.state.pwr_thermal_from_hvac =
                    (-self.hvac.state.pwr_p - self.hvac.state.pwr_i - self.hvac.state.pwr_d)
                        .max(-self.hvac.pwr_thermal_max);

                if (-self.state.pwr_thermal_from_hvac / self.hvac.state.cop) > self.hvac.pwr_aux_max
                {
                    self.hvac.state.pwr_aux = self.hvac.pwr_aux_max;
                    // correct if limit is exceeded
                    self.state.pwr_thermal_from_hvac =
                        -self.hvac.state.pwr_aux * self.hvac.state.cop;
                }
            } else {
                // HEATING MODE; cabin is colder than set point

                if self.hvac.state.pwr_i > si::Power::ZERO {
                    // reset to switch from cooling to heating
                    self.hvac.state.pwr_i = si::Power::ZERO;
                }
                self.hvac.state.pwr_i = self.hvac.state.pwr_i.max(-self.hvac.pwr_i_max);

                self.state.pwr_thermal_from_hvac =
                    (-self.hvac.state.pwr_p - self.hvac.state.pwr_i - self.hvac.state.pwr_d)
                        .min(self.hvac.pwr_thermal_max);

                // Assumes blower has negligible impact on aux load, may want to revise later
                match self.hvac.heat_source {
                    HeatSource::FuelConverter => {
                        ensure!(
                            te_fc.is_some(),
                            "{}\nExpected vehicle with [FuelConverter] with thermal plant model.",
                            format_dbg!()
                        );
                        // limit heat transfer to be substantially less than what is physically possible
                        // i.e. the engine can't drop below cabin temperature to heat the cabin
                        self.state.pwr_thermal_from_hvac = self
                            .state
                            .pwr_thermal_from_hvac
                            .min(
                                self.heat_capacitance *
                                (te_fc.unwrap() - self.state.temp)
                                    * 0.1 // so that it's substantially less
                                    / dt,
                            )
                            .max(si::Power::ZERO);
                        // TODO: think about what to do for PHEV, which needs careful consideration here
                        // HEV probably also needs careful consideration
                        // There needs to be an engine temperature (e.g. 60°C) below which the engine is forced on
                    }
                    HeatSource::ResistanceHeater => {}
                    HeatSource::HeatPump => {}
                }
            }
        }
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
    // TODO: make sure this gets updated
    temp: si::Temperature,
    /// lumped cabin temperature at previous simulation time step
    // TODO: make sure this gets updated
    temp_prev: si::Temperature,
    /// Thermal power coming to cabin from HVAC system.  Positive indicates
    /// heating, and negative indicates cooling.
    pwr_thermal_from_hvac: si::Power,
    /// Cumulative thermal energy coming to cabin from HVAC system.  Positive indicates
    /// heating, and negative indicates cooling.
    energy_thermal_from_hvac: si::Energy,
    /// Thermal power coming to cabin from ambient air.  Positive indicates
    /// heating, and negative indicates cooling.
    pwr_thermal_from_amb: si::Power,
    /// Cumulative thermal energy coming to cabin from ambient air.  Positive indicates
    /// heating, and negative indicates cooling.
    energy_thermal_from_amb: si::Energy,
}

impl Init for LumpedCabinState {}
impl SerdeAPI for LumpedCabinState {}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// HVAC system for [LumpedCabin]
pub struct HVACSystemForSCC {
    /// set point temperature
    pub te_set: si::Temperature,
    /// deadband range.  any cabin temperature within this range of
    /// `te_set` results in no HVAC power draw
    pub te_deadband: si::Temperature,
    /// HVAC proportional gain
    pub p: si::ThermalConductance,
    /// HVAC integral gain [W / K / s], resets at zero crossing events
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub i: f64,
    /// value at which [Self::i] stops accumulating
    pub pwr_i_max: si::Power,
    /// HVAC derivative gain [W / K * s]  
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub d: f64,
    /// max HVAC thermal power
    pub pwr_thermal_max: si::Power,
    /// coefficient between 0 and 1 to calculate HVAC efficiency by multiplying by
    /// coefficient of performance (COP)
    pub frac_of_ideal_cop: f64,
    // TODO: make sure this is plumbed up
    /// heat source
    #[api(skip_get, skip_set)]
    pub heat_source: HeatSource,
    /// max allowed aux load
    pub pwr_aux_max: si::Power,
    /// coefficient of performance of vapor compression cycle
    pub state: HVACSystemForSCCState,
    #[serde(default)]
    #[serde(skip_serializing_if = "HVACSystemForSCCStateHistoryVec::is_empty")]
    pub history: HVACSystemForSCCStateHistoryVec,
}
impl Init for HVACSystemForSCC {}
impl SerdeAPI for HVACSystemForSCC {}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub enum HeatSource {
    /// [FuelConverter], if applicable, provides heat for HVAC system
    FuelConverter,
    /// Resistance heater provides heat for HVAC system
    ResistanceHeater,
    /// Heat pump provides heat for HVAC system
    HeatPump,
}
impl Init for HeatSource {}
impl SerdeAPI for HeatSource {}

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
pub struct HVACSystemForSCCState {
    /// time step counter
    pub i: u32,
    /// portion of total HVAC cooling/heating (negative/positive) power due to proportional gain
    pwr_p: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to proportional gain
    energy_p: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power due to integral gain
    pwr_i: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to integral gain
    energy_i: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power due to derivative gain
    pwr_d: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to derivative gain
    energy_d: si::Energy,
    /// input power required
    pwr_in: si::Power,
    /// input cumulative energy required
    energy_in: si::Energy,
    /// coefficient of performance (i.e. efficiency) of vapor compression cycle
    cop: si::Ratio,
    /// Aux power demand from HVAC system
    pwr_aux: si::Power,
    /// Cumulative aux energy for HVAC system
    energy_aux: si::Energy,
}
impl Init for HVACSystemForSCCState {}
impl SerdeAPI for HVACSystemForSCCState {}
