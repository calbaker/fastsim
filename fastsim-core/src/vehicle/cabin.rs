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
    /// HVAC model
    pub hvac: HVACSystemForSCC,
    pub state: LumpedCabinState,
    #[serde(default)]
    #[serde(skip_serializing_if = "LumpedCabinStateHistoryVec::is_empty")]
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
    /// - `te_fc`: [FuelConverter] temperature, as appropriate for [PowertrainType]
    /// - `dt`: simulation time step size
    pub fn solve(
        &mut self,
        te_amb_air: si::Temperature,
        te_fc: Option<si::Temperature>,
        veh_state: VehicleState,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        self.state.pwr_thermal_from_hvac = self.hvac.get_pwr_thermal_from_hvac(
            te_amb_air,
            te_fc,
            self.state,
            self.heat_capacitance,
            dt,
        )?;

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
    // TODO: make sure this gets updated
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
    /// heat source
    #[api(skip_get, skip_set)]
    pub heat_source: CabinHeatSource,
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
impl HVACSystemForSCC {
    pub fn get_pwr_thermal_from_hvac(
        &mut self,
        te_amb_air: si::Temperature,
        te_fc: Option<si::Temperature>,
        cab_state: LumpedCabinState,
        cab_heat_cap: si::HeatCapacity,
        dt: si::Time,
    ) -> anyhow::Result<si::Power> {
        let pwr_from_hvac = if cab_state.temp <= self.te_set + self.te_deadband
            && cab_state.temp >= self.te_set - self.te_deadband
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_i = si::Power::ZERO; // reset to 0.0
            self.state.pwr_p = si::Power::ZERO;
            self.state.pwr_d = si::Power::ZERO;
            si::Power::ZERO
        } else {
            let te_delta_vs_set = cab_state.temp - self.te_set;
            let te_delta_vs_amb: si::Temperature = cab_state.temp - te_amb_air;

            self.state.pwr_p = -self.p * te_delta_vs_set;
            self.state.pwr_i -= self.i * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt;
            self.state.pwr_i = self.state.pwr_i.max(-self.pwr_i_max).min(self.pwr_i_max);
            self.state.pwr_d =
                -self.d * uc::J / uc::KELVIN * ((cab_state.temp - cab_state.temp_prev) / dt);

            // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
            // cop_ideal is t_h / (t_h - t_c) for heating
            // cop_ideal is t_c / (t_h - t_c) for cooling

            // divide-by-zero protection and realistic limit on COP
            let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                // cabin is cooler than ambient + threshold
                // TODO: make this `5.0` not hardcoded
                cab_state.temp / (5.0 * uc::KELVIN)
            } else {
                cab_state.temp / te_delta_vs_amb.abs()
            };
            self.state.cop = cop_ideal * self.frac_of_ideal_cop;
            assert!(self.state.cop > 0.0 * uc::R);

            if cab_state.temp > self.te_set + self.te_deadband {
                // COOLING MODE; cabin is hotter than set point

                if self.state.pwr_i > si::Power::ZERO {
                    // If `pwr_i` is greater than zero, reset to switch from heating to cooling
                    self.state.pwr_i = si::Power::ZERO;
                }
                let mut pwr_thermal_from_hvac =
                    (self.state.pwr_p + self.state.pwr_i + self.state.pwr_d)
                        .max(-self.pwr_thermal_max);

                if (-pwr_thermal_from_hvac / self.state.cop) > self.pwr_aux_max {
                    self.state.pwr_aux = self.pwr_aux_max;
                    // correct if limit is exceeded
                    pwr_thermal_from_hvac = -self.state.pwr_aux * self.state.cop;
                } else {
                    self.state.pwr_aux = pwr_thermal_from_hvac / self.state.cop;
                }
                self.state.pwr_thermal_req = si::Power::ZERO;
                pwr_thermal_from_hvac
            } else {
                // HEATING MODE; cabin is colder than set point

                if self.state.pwr_i < si::Power::ZERO {
                    // If `pwr_i` is less than zero reset to switch from cooling to heating
                    self.state.pwr_i = si::Power::ZERO;
                }
                let mut pwr_thermal_from_hvac =
                    (-self.state.pwr_p - self.state.pwr_i - self.state.pwr_d)
                        .min(self.pwr_thermal_max);

                // Assumes blower has negligible impact on aux load, may want to revise later
                match self.heat_source {
                    CabinHeatSource::FuelConverter => {
                        ensure!(
                            te_fc.is_some(),
                            "{}\nExpected vehicle with [FuelConverter] with thermal plant model.",
                            format_dbg!()
                        );
                        // limit heat transfer to be substantially less than what is physically possible
                        // i.e. the engine can't drop below cabin temperature to heat the cabin
                        pwr_thermal_from_hvac = pwr_thermal_from_hvac
                            .min(
                                cab_heat_cap *
                                (te_fc.unwrap() - cab_state.temp)
                                    * 0.1 // so that it's substantially less
                                    / dt,
                            )
                            .max(si::Power::ZERO);
                        self.state.cop = f64::NAN * uc::R;
                        self.state.pwr_thermal_req = pwr_thermal_from_hvac;
                        // Assumes aux power needed for heating is incorporated into based aux load.
                        // TODO: refine this, perhaps by making aux power
                        // proportional to heating power, to account for blower power
                        self.state.pwr_aux = si::Power::ZERO;
                        // TODO: think about what to do for PHEV, which needs careful consideration here
                        // HEV probably also needs careful consideration
                        // There needs to be an engine temperature (e.g. 60°C) below which the engine is forced on
                    }
                    CabinHeatSource::ResistanceHeater => {
                        self.state.cop = uc::R;
                        self.state.pwr_thermal_req = si::Power::ZERO;
                        self.state.pwr_aux = pwr_thermal_from_hvac;
                    }
                    CabinHeatSource::HeatPump => {
                        self.state.pwr_thermal_req = si::Power::ZERO;
                        // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
                        // cop_ideal is t_h / (t_h - t_c) for heating
                        // cop_ideal is t_c / (t_h - t_c) for cooling

                        // divide-by-zero protection and realistic limit on COP
                        // TODO: make sure this is right for heating!
                        let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                            // cabin is cooler than ambient + threshold
                            // TODO: make this `5.0` not hardcoded
                            cab_state.temp / (5.0 * uc::KELVIN)
                        } else {
                            cab_state.temp / te_delta_vs_amb.abs()
                        };
                        self.state.cop = cop_ideal * self.frac_of_ideal_cop;
                        assert!(self.state.cop > 0.0 * uc::R);
                        self.state.pwr_thermal_req = si::Power::ZERO;
                        if (pwr_thermal_from_hvac / self.state.cop) > self.pwr_aux_max {
                            self.state.pwr_aux = self.pwr_aux_max;
                            // correct if limit is exceeded
                            pwr_thermal_from_hvac = -self.state.pwr_aux * self.state.cop;
                        } else {
                            self.state.pwr_aux = pwr_thermal_from_hvac / self.state.cop;
                        }
                    }
                }
                pwr_thermal_from_hvac
            }
        };
        Ok(pwr_from_hvac)
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
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

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
pub struct HVACSystemForSCCState {
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
    /// Aux power demand from HVAC system
    pub pwr_aux: si::Power,
    /// Cumulative aux energy for HVAC system
    pub energy_aux: si::Energy,
    /// Thermal power demand by HVAC system from thermal component (e.g. [FuelConverter])
    pub pwr_thermal_req: si::Power,
    /// Cumulative energy demand by HVAC system from thermal component (e.g. [FuelConverter])
    pub energy_thermal_req: si::Energy,
}
impl Init for HVACSystemForSCCState {}
impl SerdeAPI for HVACSystemForSCCState {}
