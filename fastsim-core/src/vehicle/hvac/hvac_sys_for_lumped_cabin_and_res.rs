use super::*;

#[fastsim_api(
    #[staticmethod]
    #[pyo3(name = "default")]
    fn default_py() -> Self {
        Default::default()
    }
)]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// HVAC system for [LumpedCabin] and [ReversibleEnergyStorage::thrml]
pub struct HVACSystemForLumpedCabinAndRES {
    /// set point temperature
    pub te_set: si::Temperature,
    /// Deadband range.  Any cabin temperature within this range of `te_set`
    /// results in no HVAC power draw
    pub te_deadband: si::TemperatureInterval,
    /// HVAC proportional gain for cabin
    pub p_cabin: si::ThermalConductance,
    /// HVAC integral gain [W / K / s] for cabin, resets at zero crossing events  
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub i_cabin: f64,
    /// value at which state.i stops accumulating for cabin
    pub pwr_i_max_cabin: si::Power,
    /// HVAC derivative gain [W / K * s] for cabin
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub d_cabin: f64,
    /// max HVAC thermal power
    /// HVAC proportional gain for [ReversibleEnergyStorage]
    pub p_res: si::ThermalConductance,
    /// HVAC integral gain [W / K / s] for [ReversibleEnergyStorage], resets at zero crossing events  
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub i_res: f64,
    /// value at which state.i stops accumulating for [ReversibleEnergyStorage]
    pub pwr_i_max_res: si::Power,
    /// HVAC derivative gain [W / K * s] for [ReversibleEnergyStorage]
    /// NOTE: `uom` crate does not have this unit, but it may be possible to make a custom unit for this
    pub d_res: f64,
    /// max HVAC thermal power
    pub pwr_thrml_max: si::Power,
    /// coefficient between 0 and 1 to calculate HVAC efficiency by multiplying by
    /// coefficient of performance (COP)
    pub frac_of_ideal_cop: f64,
    /// cabin heat source
    pub cabin_heat_source: CabinHeatSource,
    /// res heat source
    pub res_heat_source: RESHeatSource,
    /// res cooling source
    pub res_cooling_source: RESCoolingSource,
    /// max allowed aux load
    pub pwr_aux_for_hvac_max: si::Power,
    /// coefficient of performance of vapor compression cycle
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: HVACSystemForLumpedCabinAndRESState,
    #[serde(
        default,
        skip_serializing_if = "HVACSystemForLumpedCabinAndRESStateHistoryVec::is_empty"
    )]
    pub history: HVACSystemForLumpedCabinAndRESStateHistoryVec,
}
impl Default for HVACSystemForLumpedCabinAndRES {
    fn default() -> Self {
        Self {
            te_set: *TE_STD_AIR,
            te_deadband: 1.5 * uc::KELVIN_INT,
            p_cabin: Default::default(),
            i_cabin: Default::default(),
            d_cabin: Default::default(),
            pwr_i_max_cabin: 5. * uc::KW,
            p_res: Default::default(),
            i_res: Default::default(),
            d_res: Default::default(),
            pwr_i_max_res: 5. * uc::KW,
            pwr_thrml_max: 10. * uc::KW,
            frac_of_ideal_cop: 0.15,
            cabin_heat_source: CabinHeatSource::ResistanceHeater,
            res_heat_source: RESHeatSource::ResistanceHeater,
            res_cooling_source: RESCoolingSource::HVAC,
            pwr_aux_for_hvac_max: uc::KW * 5.,
            state: Default::default(),
            history: Default::default(),
        }
    }
}
impl SetCumulative for HVACSystemForLumpedCabinAndRES {
    fn set_cumulative(&mut self, dt: si::Time) {
        self.state.set_cumulative(dt);
    }
}
impl Init for HVACSystemForLumpedCabinAndRES {}
impl SerdeAPI for HVACSystemForLumpedCabinAndRES {}
impl HVACSystemForLumpedCabinAndRES {
    /// # Arguments
    /// - `te_amb_air`: ambient air temperature
    /// - `te_fc`: [FuelConverter] temperature, if equipped
    /// - `cab_state`: [LumpedCabinState]
    /// - `cab_heat_cap`: [LumpedCabinState] heat capacity
    /// - `res_temp`: [ReversibleEnergyStorage] temperatures at current and previous time step
    /// - `dt`: simulation time step size
    ///
    /// # Returns
    /// - `pwr_thrml_hvac_to_cabin`: thermal power flowing from [Vehicle::hvac] system to cabin  
    /// - `pwr_thrml_fc_to_cabin`: thermal power flowing from [FuelConverter] to cabin  
    /// - `pwr_thrml_hvac_to_res`: thermal power flowing from [Vehicle::hvac] system to
    ///     [ReversibleEnergyStorage] `thrml` system  
    #[allow(clippy::too_many_arguments)] // the order is reasonably protected by typing
    pub fn solve(
        &mut self,
        te_amb_air: si::Temperature,
        te_fc: Option<si::Temperature>,
        cab_state: LumpedCabinState,
        cab_heat_cap: si::HeatCapacity,
        res_temps: (si::Temperature, si::Temperature),
        dt: si::Time,
    ) -> anyhow::Result<(si::Power, si::Power, si::Power)> {
        let (res_temp, res_temp_prev) = res_temps;
        ensure!(!res_temp.is_nan(), format_dbg!(res_temp));
        ensure!(!res_temp_prev.is_nan(), format_dbg!(res_temp_prev));
        let mut pwr_thrml_hvac_to_cabin = self
            .solve_for_cabin(te_fc, cab_state, cab_heat_cap, dt)
            .with_context(|| format_dbg!())?;
        let mut pwr_thrml_hvac_to_res: si::Power = self
            .solve_for_res(res_temp, res_temp_prev, dt)
            .with_context(|| format_dbg!())?;
        let (cop_ideal, te_ref) = if pwr_thrml_hvac_to_res + pwr_thrml_hvac_to_cabin
            > si::Power::ZERO
        {
            // heating mode
            // TODO: account for cabin and battery heat sources in COP calculation!!!!

            let (te_ref, te_delta_vs_amb) = if pwr_thrml_hvac_to_res > si::Power::ZERO {
                // both powers are positive -- i.e. both are in heating mode

                let te_ref: si::Temperature = if cab_state.temperature > res_temp {
                    // cabin is hotter
                    cab_state.temperature
                } else {
                    // battery is hotter
                    res_temp
                };
                (
                    te_ref,
                    (te_ref.get::<si::degree_celsius>() - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            } else if pwr_thrml_hvac_to_res >= si::Power::ZERO {
                // `pwr_thrml_hvac_to_res` dominates need for heating
                (
                    res_temp,
                    (res_temp.get::<si::degree_celsius>() - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            } else {
                // `pwr_thrml_hvac_to_res` dominates need for heating
                (
                    cab_state.temperature,
                    (cab_state.temperature.get::<si::degree_celsius>()
                        - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            };

            // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
            // cop_ideal is t_h / (t_h - t_c) for heating
            // cop_ideal is t_c / (t_h - t_c) for cooling

            // divide-by-zero protection and realistic limit on COP
            // TODO: make sure this is consistent with above commented equation for heating!
            if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN_INT {
                // cabin is cooler than ambient + threshold
                // TODO: make this `5.0` not hardcoded
                (te_ref / (5.0 * uc::KELVIN), te_ref)
            } else {
                (te_ref / te_delta_vs_amb.abs(), te_ref)
            }
        } else if pwr_thrml_hvac_to_res + pwr_thrml_hvac_to_cabin < si::Power::ZERO {
            // cooling mode
            // TODO: account for battery cooling source in COP calculation!!!!

            let (te_ref, te_delta_vs_amb) = if pwr_thrml_hvac_to_res < si::Power::ZERO {
                // both powers are negative -- i.e. both are in cooling mode

                let te_ref: si::Temperature = if cab_state.temperature < res_temp {
                    // cabin is colder
                    cab_state.temperature
                } else {
                    // battery is colder
                    res_temp
                };
                (
                    te_ref,
                    (te_ref.get::<si::degree_celsius>() - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            } else if pwr_thrml_hvac_to_res >= si::Power::ZERO {
                // `pwr_thrml_hvac_to_cabin` dominates need for cooling
                (
                    cab_state.temperature,
                    (cab_state.temperature.get::<si::kelvin_abs>()
                        - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            } else {
                // `pwr_thrml_hvac_to_res` dominates need for cooling
                (
                    res_temp,
                    (res_temp.get::<si::degree_celsius>() - te_amb_air.get::<si::degree_celsius>())
                        * uc::KELVIN_INT,
                )
            };

            // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
            // cop_ideal is t_h / (t_h - t_c) for heating
            // cop_ideal is t_c / (t_h - t_c) for cooling

            // divide-by-zero protection and realistic limit on COP
            if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN_INT {
                // cooling-dominating component is cooler than ambient + threshold
                // TODO: make this `5.0` not hardcoded
                (te_ref / (5.0 * uc::KELVIN), te_ref)
            } else {
                (te_ref / te_delta_vs_amb.abs(), te_ref)
            }
        } else {
            (si::Ratio::ZERO, f64::NAN * uc::KELVIN)
        };
        self.state.cop = cop_ideal * self.frac_of_ideal_cop;
        ensure!(
            self.state.cop >= 0.0 * uc::R,
            "{}\n{}",
            format_dbg!(cop_ideal),
            format_dbg!(te_ref)
        );

        let mut pwr_thrml_fc_to_cabin = si::Power::ZERO;
        self.state.pwr_aux_for_hvac = if pwr_thrml_hvac_to_cabin > si::Power::ZERO {
            match self.cabin_heat_source {
                CabinHeatSource::FuelConverter => {
                    pwr_thrml_fc_to_cabin = pwr_thrml_hvac_to_cabin;
                    // NOTE: should make this scale with power demand
                    si::Power::ZERO
                }
                CabinHeatSource::ResistanceHeater => pwr_thrml_hvac_to_cabin,
                CabinHeatSource::HeatPump => pwr_thrml_hvac_to_cabin * self.state.cop,
            }
        } else {
            -pwr_thrml_hvac_to_cabin * self.state.cop
        } + if pwr_thrml_hvac_to_res > si::Power::ZERO {
            match self.res_heat_source {
                RESHeatSource::ResistanceHeater => pwr_thrml_hvac_to_res,
                RESHeatSource::HeatPump => pwr_thrml_hvac_to_res * self.state.cop,
                RESHeatSource::None => {
                    pwr_thrml_hvac_to_res = si::Power::ZERO;
                    si::Power::ZERO
                }
            }
        } else {
            match self.res_cooling_source {
                RESCoolingSource::HVAC => -pwr_thrml_hvac_to_res * self.state.cop,
                RESCoolingSource::None => {
                    pwr_thrml_hvac_to_res = si::Power::ZERO;
                    si::Power::ZERO
                }
            }
        };

        self.state.pwr_aux_for_hvac = if self.state.pwr_aux_for_hvac > self.pwr_aux_for_hvac_max {
            pwr_thrml_hvac_to_res =
                self.pwr_aux_for_hvac_max * self.state.cop * pwr_thrml_hvac_to_res
                    / (pwr_thrml_hvac_to_res + pwr_thrml_hvac_to_cabin);
            pwr_thrml_hvac_to_cabin =
                self.pwr_aux_for_hvac_max * self.state.cop * pwr_thrml_hvac_to_cabin
                    / (pwr_thrml_hvac_to_res + pwr_thrml_hvac_to_cabin);
            if pwr_thrml_hvac_to_cabin > si::Power::ZERO
                && self.cabin_heat_source.is_fuel_converter()
            {
                pwr_thrml_fc_to_cabin = pwr_thrml_hvac_to_cabin;
            }
            self.pwr_aux_for_hvac_max
        } else {
            self.state.pwr_aux_for_hvac
        };

        self.state.pwr_thrml_hvac_to_cabin = pwr_thrml_hvac_to_cabin;
        self.state.pwr_thrml_fc_to_cabin = pwr_thrml_fc_to_cabin;
        self.state.pwr_thrml_hvac_to_res = pwr_thrml_hvac_to_res;

        Ok((
            self.state.pwr_thrml_hvac_to_cabin,
            self.state.pwr_thrml_fc_to_cabin,
            self.state.pwr_thrml_hvac_to_res,
        ))
    }

    fn solve_for_cabin(
        &mut self,
        te_fc: Option<si::Temperature>,
        cab_state: LumpedCabinState,
        cab_heat_cap: si::HeatCapacity,
        dt: si::Time,
    ) -> anyhow::Result<si::Power> {
        let pwr_thrml_hvac_to_cabin = if cab_state.temperature <= (self.te_set + self.te_deadband)
            && cab_state.temperature >= (self.te_set - self.te_deadband)
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_i = si::Power::ZERO; // reset to 0.0
            self.state.pwr_p = si::Power::ZERO;
            self.state.pwr_d = si::Power::ZERO;
            si::Power::ZERO
        } else {
            // outside deadband
            let te_delta_vs_set = (cab_state.temperature.get::<si::degree_celsius>()
                - self.te_set.get::<si::degree_celsius>())
                * uc::KELVIN_INT;

            self.state.pwr_p = -self.p_cabin * te_delta_vs_set;
            self.state.pwr_i -= self.i_cabin * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt;
            self.state.pwr_i = self
                .state
                .pwr_i
                .max(-self.pwr_i_max_cabin)
                .min(self.pwr_i_max_cabin);
            self.state.pwr_d = -self.d_cabin * uc::J / uc::KELVIN
                * ((cab_state.temperature.get::<si::degree_celsius>()
                    - cab_state.temp_prev.get::<si::degree_celsius>())
                    * uc::KELVIN_INT
                    / dt);

            let pwr_thrml_hvac_to_cabin: si::Power =
                if cab_state.temperature > self.te_set + self.te_deadband {
                    // COOLING MODE; cabin is hotter than set point

                    if self.state.pwr_i > si::Power::ZERO {
                        // If `pwr_i` is greater than zero, reset to switch from heating to cooling
                        self.state.pwr_i = si::Power::ZERO;
                    }
                    let mut pwr_thrml_hvac_to_cab =
                        (self.state.pwr_p + self.state.pwr_i + self.state.pwr_d)
                            .max(-self.pwr_thrml_max);

                    if (-pwr_thrml_hvac_to_cab / self.state.cop) > self.pwr_aux_for_hvac_max {
                        self.state.pwr_aux_for_hvac = self.pwr_aux_for_hvac_max;
                        // correct if limit is exceeded
                        pwr_thrml_hvac_to_cab = -self.state.pwr_aux_for_hvac * self.state.cop;
                    } else {
                        self.state.pwr_aux_for_hvac = pwr_thrml_hvac_to_cab / self.state.cop;
                    }
                    pwr_thrml_hvac_to_cab
                } else {
                    // HEATING MODE; cabin is colder than set point

                    if self.state.pwr_i < si::Power::ZERO {
                        // If `pwr_i` is less than zero reset to switch from cooling to heating
                        self.state.pwr_i = si::Power::ZERO;
                    }
                    let mut pwr_thrml_hvac_to_cabin: si::Power =
                        (-self.state.pwr_p - self.state.pwr_i - self.state.pwr_d)
                            .min(self.pwr_thrml_max);

                    // Assumes blower has negligible impact on aux load, may want to revise later
                    self.handle_cabin_heat_source(
                        te_fc,
                        &mut pwr_thrml_hvac_to_cabin,
                        cab_heat_cap,
                        cab_state,
                        dt,
                    )
                    .with_context(|| format_dbg!())?;
                    pwr_thrml_hvac_to_cabin
                };
            pwr_thrml_hvac_to_cabin
        };
        Ok(pwr_thrml_hvac_to_cabin)
    }

    fn solve_for_res(
        &mut self,
        // reversible energy storage temp
        res_temp: si::Temperature,
        // reversible energy storage temp at previous time step
        res_temp_prev: si::Temperature,
        dt: si::Time,
    ) -> anyhow::Result<si::Power> {
        let pwr_thrml_hvac_to_res = if res_temp <= self.te_set + self.te_deadband
            && res_temp >= self.te_set - self.te_deadband
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_i_res = si::Power::ZERO; // reset to 0.0
            self.state.pwr_p_res = si::Power::ZERO;
            self.state.pwr_d_res = si::Power::ZERO;
            si::Power::ZERO
        } else {
            // outside deadband
            let te_delta_vs_set = (res_temp.get::<si::degree_celsius>()
                - self.te_set.get::<si::degree_celsius>())
                * uc::KELVIN_INT;
            self.state.pwr_p_res = -self.p_res * te_delta_vs_set;
            self.state.pwr_i_res -= self.i_res * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt;
            self.state.pwr_i_res = self
                .state
                .pwr_i_res
                .max(-self.pwr_i_max_res)
                .min(self.pwr_i_max_res);
            self.state.pwr_d_res = -self.d_res * uc::J / uc::KELVIN
                * ((res_temp.get::<si::degree_celsius>()
                    - res_temp_prev.get::<si::degree_celsius>())
                    * uc::KELVIN_INT
                    / dt);

            let pwr_thrml_hvac_to_res: si::Power = if res_temp > self.te_set + self.te_deadband {
                // COOLING MODE; Reversible Energy Storage is hotter than set point

                if self.state.pwr_i_res > si::Power::ZERO {
                    // If `pwr_i_res` is greater than zero, reset to switch from heating to cooling
                    self.state.pwr_i_res = si::Power::ZERO;
                }
                let mut pwr_thrml_hvac_to_res =
                    (self.state.pwr_p_res + self.state.pwr_i_res + self.state.pwr_d_res)
                        .max(-self.pwr_thrml_max);

                if (-pwr_thrml_hvac_to_res / self.state.cop) > self.pwr_aux_for_hvac_max {
                    self.state.pwr_aux_for_hvac = self.pwr_aux_for_hvac_max;
                    // correct if limit is exceeded
                    pwr_thrml_hvac_to_res = -self.state.pwr_aux_for_hvac * self.state.cop;
                } else {
                    self.state.pwr_aux_for_hvac = pwr_thrml_hvac_to_res / self.state.cop;
                }
                pwr_thrml_hvac_to_res
            } else {
                // HEATING MODE; Reversible Energy Storage is colder than set point

                if self.state.pwr_i_res < si::Power::ZERO {
                    // If `pwr_i_res` is less than zero reset to switch from cooling to heating
                    self.state.pwr_i_res = si::Power::ZERO;
                }
                #[allow(clippy::let_and_return)] // for readability
                let pwr_thrml_hvac_to_res =
                    (-self.state.pwr_p_res - self.state.pwr_i_res - self.state.pwr_d_res)
                        .min(self.pwr_thrml_max);
                pwr_thrml_hvac_to_res
            };
            pwr_thrml_hvac_to_res
        };
        Ok(pwr_thrml_hvac_to_res)
    }

    fn handle_cabin_heat_source(
        &mut self,
        te_fc: Option<si::Temperature>,
        pwr_thrml_hvac_to_cabin: &mut si::Power,
        cab_heat_cap: si::HeatCapacity,
        cab_state: LumpedCabinState,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        match self.cabin_heat_source {
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
                    (te_fc.unwrap().get::<si::degree_celsius>() - cab_state.temperature.get::<si::degree_celsius>()) * uc::KELVIN_INT
                        * 0.1 // so that it's substantially less
                        / dt,
                    )
                    .max(si::Power::ZERO);
            }
            CabinHeatSource::ResistanceHeater => {
                *pwr_thrml_hvac_to_cabin = si::Power::ZERO;
            }
            CabinHeatSource::HeatPump => {
                *pwr_thrml_hvac_to_cabin = si::Power::ZERO;
            }
        };
        Ok(())
    }
}

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
#[serde(default)]
pub struct HVACSystemForLumpedCabinAndRESState {
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
    /// portion of total HVAC cooling/heating (negative/positive) power to [ReversibleEnergyStorage::thrml] due to proportional gain
    pub pwr_p_res: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy to [ReversibleEnergyStorage::thrml] due to proportional gain
    pub energy_p_res: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power to [ReversibleEnergyStorage::thrml] due to integral gain
    pub pwr_i_res: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy to [ReversibleEnergyStorage::thrml] due to integral gain
    pub energy_i_res: si::Energy,
    /// portion of total HVAC cooling/heating (negative/positive) power to [ReversibleEnergyStorage::thrml] due to derivative gain
    pub pwr_d_res: si::Power,
    /// portion of total HVAC cooling/heating (negative/positive) cumulative energy to [ReversibleEnergyStorage::thrml] due to derivative gain
    pub energy_d_res: si::Energy,
    /// coefficient of performance (i.e. efficiency) of vapor compression cycle
    pub cop: si::Ratio,
    /// Au power demand from [Vehicle::hvac] system
    pub pwr_aux_for_hvac: si::Power,
    /// Cumulative aux energy for HVAC system
    pub energy_aux_for_hvac: si::Energy,
    /// Thermal power from HVAC system to cabin, positive is heating the cabin
    pub pwr_thrml_hvac_to_cabin: si::Power,
    /// Cumulative thermal energy from HVAC system to cabin, positive is heating the cabin
    pub energy_thrml_hvac_to_cabin: si::Energy,
    /// Thermal power from [FuelConverter] to [Cabin]
    pub pwr_thrml_fc_to_cabin: si::Power,
    /// Cumulative thermal energy from [FuelConverter] to [Cabin]
    pub energy_thrml_fc_to_cabin: si::Energy,
    /// Thermal power from HVAC to [ReversibleEnergyStorage]
    pub pwr_thrml_hvac_to_res: si::Power,
    /// Cumulative thermal energy from HVAC to [ReversibleEnergyStorage]
    pub energy_thrml_hvac_to_res: si::Energy,
}
impl Init for HVACSystemForLumpedCabinAndRESState {}
impl SerdeAPI for HVACSystemForLumpedCabinAndRESState {}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, IsVariant, From, TryInto)]
/// Heat source for [RESLumpedThermal]
pub enum RESHeatSource {
    /// Resistance heater provides heat for HVAC system
    ResistanceHeater,
    /// Heat pump provides heat for HVAC system
    HeatPump,
    /// The battery is not actively heated
    None,
}
impl Init for RESHeatSource {}
impl SerdeAPI for RESHeatSource {}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, IsVariant, From, TryInto)]
/// Cooling source for [RESLumpedThermal]
pub enum RESCoolingSource {
    /// Vapor compression system used for cabin HVAC also cools [RESLumpedThermal]
    HVAC,
    /// [RESLumpedThermal] is not actively cooled
    None,
}
impl Init for RESCoolingSource {}
impl SerdeAPI for RESCoolingSource {}
