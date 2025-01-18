use super::*;

#[fastsim_api(
    #[staticmethod]
    #[pyo3(name = "default")]
    fn default_py() -> Self {
        Default::default()
    }
)]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// HVAC system for [LumpedCabin]
pub struct HVACSystemForLumpedCabin {
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
    /// value at which state.i stops accumulating
    pub pwr_i_max: si::Power,
    /// HVAC derivative gain [W / K * s]  
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub d: f64,
    /// max HVAC thermal power
    pub pwr_thermal_max: si::Power,
    /// coefficient between 0 and 1 to calculate HVAC efficiency by multiplying by
    /// coefficient of performance (COP)
    pub frac_of_ideal_cop: f64,
    /// heat source
    pub heat_source: CabinHeatSource,
    /// max allowed aux load for HVAC
    pub pwr_aux_for_hvac_max: si::Power,
    /// coefficient of performance of vapor compression cycle
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: HVACSystemForLumpedCabinState,
    #[serde(
        default,
        skip_serializing_if = "HVACSystemForLumpedCabinStateHistoryVec::is_empty"
    )]
    pub history: HVACSystemForLumpedCabinStateHistoryVec,
}
impl Default for HVACSystemForLumpedCabin {
    fn default() -> Self {
        Self {
            te_set: *TE_STD_AIR,
            te_deadband: 1.5 * uc::KELVIN,
            p: Default::default(),
            i: Default::default(),
            d: Default::default(),
            pwr_i_max: 5. * uc::KW,
            pwr_thermal_max: 10. * uc::KW,
            frac_of_ideal_cop: 0.15,
            heat_source: CabinHeatSource::ResistanceHeater,
            pwr_aux_for_hvac_max: uc::KW * 5.,
            state: Default::default(),
            history: Default::default(),
        }
    }
}
impl Init for HVACSystemForLumpedCabin {}
impl SerdeAPI for HVACSystemForLumpedCabin {}
impl HVACSystemForLumpedCabin {
    /// # Arguments
    /// - `te_amb_air`: ambient air temperature
    /// - `te_fc`: [FuelConverter] temperature, if equipped
    /// - `cab_state`: [LumpedCabinState]
    /// - `cab_heat_cap`: cabin heat capacitance
    /// - `dt`: simulation time step size
    /// # Returns
    /// - `pwr_thrml_hvac_to_cabin`: thermal power flow from [Vehicle::hvac] system to cabin
    /// - `pwr_thrml_fc_to_cabin`: thermal power flow from [FuelConverter] to cabin via HVAC system
    pub fn solve(
        &mut self,
        te_amb_air: si::Temperature,
        te_fc: Option<si::Temperature>,
        cab_state: LumpedCabinState,
        cab_heat_cap: si::HeatCapacity,
        dt: si::Time,
    ) -> anyhow::Result<(si::Power, si::Power)> {
        let (pwr_thrml_hvac_to_cabin, pwr_thrml_fc_to_cabin, cop) = if cab_state.temperature
            <= self.te_set + self.te_deadband
            && cab_state.temperature >= self.te_set - self.te_deadband
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_i = si::Power::ZERO; // reset to 0.0
            self.state.pwr_p = si::Power::ZERO;
            self.state.pwr_d = si::Power::ZERO;
            (si::Power::ZERO, si::Power::ZERO, f64::NAN * uc::R)
        } else {
            // outside deadband
            let te_delta_vs_set = cab_state.temperature - self.te_set;
            let te_delta_vs_amb: si::Temperature = cab_state.temperature - te_amb_air;

            self.state.pwr_p = -self.p * te_delta_vs_set;
            self.state.pwr_i -= self.i * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt;
            self.state.pwr_i = self.state.pwr_i.max(-self.pwr_i_max).min(self.pwr_i_max);
            self.state.pwr_d =
                -self.d * uc::J / uc::KELVIN * ((cab_state.temperature - cab_state.temp_prev) / dt);

            let (pwr_thrml_hvac_to_cabin, pwr_thrml_fc_to_cabin, cop) =
                if cab_state.temperature > self.te_set + self.te_deadband {
                    // COOLING MODE; cabin is hotter than set point

                    // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
                    // cop_ideal is t_h / (t_h - t_c) for heating
                    // cop_ideal is t_c / (t_h - t_c) for cooling

                    // divide-by-zero protection and realistic limit on COP
                    let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                        // cabin is cooler than ambient + threshold
                        // TODO: make this `5.0` not hardcoded
                        cab_state.temperature / (5.0 * uc::KELVIN)
                    } else {
                        cab_state.temperature / te_delta_vs_amb.abs()
                    };
                    let cop = cop_ideal * self.frac_of_ideal_cop;
                    assert!(cop > 0.0 * uc::R);

                    if self.state.pwr_i > si::Power::ZERO {
                        // If `pwr_i` is greater than zero, reset to switch from heating to cooling
                        self.state.pwr_i = si::Power::ZERO;
                    }
                    let mut pwr_thrml_hvac_to_cab =
                        (self.state.pwr_p + self.state.pwr_i + self.state.pwr_d)
                            .max(-self.pwr_thermal_max);

                    if (-pwr_thrml_hvac_to_cab / self.state.cop) > self.pwr_aux_for_hvac_max {
                        self.state.pwr_aux_for_hvac = self.pwr_aux_for_hvac_max;
                        // correct if limit is exceeded
                        pwr_thrml_hvac_to_cab = -self.state.pwr_aux_for_hvac * self.state.cop;
                    } else {
                        self.state.pwr_aux_for_hvac = -pwr_thrml_hvac_to_cab / self.state.cop;
                    }
                    let pwr_thrml_fc_to_cabin = si::Power::ZERO;
                    (pwr_thrml_hvac_to_cab, pwr_thrml_fc_to_cabin, cop)
                } else {
                    // HEATING MODE; cabin is colder than set point

                    if self.state.pwr_i < si::Power::ZERO {
                        // If `pwr_i` is less than zero reset to switch from cooling to heating
                        self.state.pwr_i = si::Power::ZERO;
                    }
                    let mut pwr_thrml_hvac_to_cabin =
                        (-self.state.pwr_p - self.state.pwr_i - self.state.pwr_d)
                            .min(self.pwr_thermal_max);

                    // Assumes blower has negligible impact on aux load, may want to revise later
                    let (pwr_thrml_fc_to_cabin, cop) = self
                        .handle_heat_source(
                            te_fc,
                            te_delta_vs_amb,
                            &mut pwr_thrml_hvac_to_cabin,
                            cab_heat_cap,
                            cab_state,
                            dt,
                        )
                        .with_context(|| format_dbg!())?;
                    (pwr_thrml_hvac_to_cabin, pwr_thrml_fc_to_cabin, cop)
                };
            (pwr_thrml_hvac_to_cabin, pwr_thrml_fc_to_cabin, cop)
        };
        self.state.cop = cop;
        Ok((pwr_thrml_hvac_to_cabin, pwr_thrml_fc_to_cabin))
    }

    fn handle_heat_source(
        &mut self,
        te_fc: Option<si::Temperature>,
        te_delta_vs_amb: si::Temperature,
        pwr_thrml_hvac_to_cabin: &mut si::Power,
        cab_heat_cap: si::HeatCapacity,
        cab_state: LumpedCabinState,
        dt: si::Time,
    ) -> anyhow::Result<(si::Power, si::Ratio)> {
        let (pwr_thrml_fc_to_cabin, cop) = match self.heat_source {
            CabinHeatSource::FuelConverter => {
                ensure!(
                    te_fc.is_some(),
                    "{}\nExpected vehicle with [FuelConverter] with thermal plant model.",
                    format_dbg!()
                );
                // limit heat transfer to be substantially less than what is physically possible
                // i.e. the engine can't drop below cabin temperature to heat the cabin
                *pwr_thrml_hvac_to_cabin = pwr_thrml_hvac_to_cabin
                    .min(
                        cab_heat_cap *
                    (te_fc.unwrap() - cab_state.temperature)
                        * 0.1 // so that it's substantially less
                        / dt,
                    )
                    .max(si::Power::ZERO);
                let cop = f64::NAN * uc::R;
                let pwr_thrml_fc_to_cabin = *pwr_thrml_hvac_to_cabin;
                // Assumes aux power needed for heating is incorporated into based aux load.
                // TODO: refine this, perhaps by making aux power
                // proportional to heating power, to account for blower power
                self.state.pwr_aux_for_hvac = si::Power::ZERO;
                (pwr_thrml_fc_to_cabin, cop)
            }
            CabinHeatSource::ResistanceHeater => {
                let cop = uc::R;
                self.state.pwr_aux_for_hvac = *pwr_thrml_hvac_to_cabin; // COP is 1 so does not matter
                #[allow(clippy::let_and_return)] // for readability
                let pwr_thrml_fc_to_cabin = si::Power::ZERO;
                (pwr_thrml_fc_to_cabin, cop)
            }
            CabinHeatSource::HeatPump => {
                // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
                // cop_ideal is t_h / (t_h - t_c) for heating
                // cop_ideal is t_c / (t_h - t_c) for cooling

                // divide-by-zero protection and realistic limit on COP
                // TODO: make sure this is consist with above commented equation for heating!
                let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                    // cabin is cooler than ambient + threshold
                    // TODO: make this `5.0` not hardcoded
                    cab_state.temperature / (5.0 * uc::KELVIN)
                } else {
                    cab_state.temperature / te_delta_vs_amb.abs()
                };
                let cop = cop_ideal * self.frac_of_ideal_cop;
                assert!(cop > 0.0 * uc::R);
                if (*pwr_thrml_hvac_to_cabin / self.state.cop) > self.pwr_aux_for_hvac_max {
                    self.state.pwr_aux_for_hvac = self.pwr_aux_for_hvac_max;
                    // correct if limit is exceeded
                    *pwr_thrml_hvac_to_cabin = -self.state.pwr_aux_for_hvac * self.state.cop;
                } else {
                    self.state.pwr_aux_for_hvac = *pwr_thrml_hvac_to_cabin / self.state.cop;
                }
                #[allow(clippy::let_and_return)] // for readability
                let pwr_thrml_fc_to_cabin = si::Power::ZERO;
                (pwr_thrml_fc_to_cabin, cop)
            }
        };
        Ok((pwr_thrml_fc_to_cabin, cop))
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, IsVariant, From, TryInto)]
pub enum CabinHeatSource {
    /// [FuelConverter], if applicable, provides heat for HVAC system
    FuelConverter,
    /// Resistance heater provides heat for HVAC system
    ResistanceHeater,
    /// Heat pump provides heat for HVAC system
    HeatPump,
}
impl Init for CabinHeatSource {}
impl SerdeAPI for CabinHeatSource {}

#[fastsim_api(
    #[pyo3(name = "default")]
    #[staticmethod]
    fn default_py() -> Self {
        Self::default()
    }
)]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
#[serde(default)]
pub struct HVACSystemForLumpedCabinState {
    /// time step counter
    pub i: u32,
    /// portion of total HVAC cooling/heating (negative/positive) power due to proportional gain
    pub pwr_p: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to proportional gain
    pub energy_p: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power due to integral gain
    pub pwr_i: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to integral gain
    pub energy_i: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power due to derivative gain
    pub pwr_d: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy due to derivative gain
    pub energy_d: si::Energy,
    /// coefficient of performance (i.e. efficiency) of vapor compression cycle
    pub cop: si::Ratio,
    /// Aux power demand from [Vehicle::hvac] system
    pub pwr_aux_for_hvac: si::Power,
    /// Cumulative aux energy for HVAC system
    pub energy_aux: si::Energy,
    /// Cumulative energy demand by HVAC system from thermal component (e.g. [FuelConverter])
    pub energy_thermal_req: si::Energy,
}
impl Init for HVACSystemForLumpedCabinState {}
impl SerdeAPI for HVACSystemForLumpedCabinState {}
