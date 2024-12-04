use super::*;
use crate::prelude::*;
use serde::Deserializer;
use std::f64::consts::PI;

// TODO: think about how to incorporate life modeling for Fuel Cells and other tech

#[fastsim_api(
    // // optional, custom, struct-specific pymethods
    // #[getter("eff_max")]
    // fn get_eff_max_py(&self) -> f64 {
    //     self.get_eff_max()
    // }

    // #[setter("__eff_max")]
    // fn set_eff_max_py(&mut self, eff_max: f64) -> PyResult<()> {
    //     self.set_eff_max(eff_max).map_err(PyValueError::new_err)
    // }

    // #[getter("eff_min")]
    // fn get_eff_min_py(&self) -> f64 {
    //     self.get_eff_min()
    // }

    // #[getter("eff_range")]
    // fn get_eff_range_py(&self) -> f64 {
    //     self.get_eff_range()
    // }

    // #[setter("__eff_range")]
    // fn set_eff_range_py(&mut self, eff_range: f64) -> PyResult<()> {
    //     self.set_eff_range(eff_range).map_err(PyValueError::new_err)
    // }

    // TODO: handle `side_effects` and uncomment
    // #[setter("__mass_kg")]
    // fn set_mass_py(&mut self, mass_kg: Option<f64>) -> anyhow::Result<()> {
    //     self.set_mass(mass_kg.map(|m| m * uc::KG))?;
    //     Ok(())
    // }

    #[getter("mass_kg")]
    fn get_mass_py(&self) -> PyResult<Option<f64>> {
        Ok(self.mass()?.map(|m| m.get::<si::kilogram>()))
    }

    #[getter]
    fn get_specific_pwr_kw_per_kg(&self) -> Option<f64> {
        self.specific_pwr.map(|x| x.get::<si::kilowatt_per_kilogram>())
    }
)]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// Struct for modeling Fuel Converter (e.g. engine, fuel cell.)
pub struct FuelConverter {
    /// [Self] Thermal plant, including thermal management controls
    #[serde(default, skip_serializing_if = "FuelConverterThermalOption::is_none")]
    #[api(skip_get, skip_set)]
    pub thrml: FuelConverterThermalOption,
    /// [Self] mass
    #[serde(default)]
    #[api(skip_get, skip_set)]
    pub(in super::super) mass: Option<si::Mass>,
    /// FuelConverter specific power
    #[api(skip_get, skip_set)]
    pub(in super::super) specific_pwr: Option<si::SpecificPower>,
    /// max rated brake output power
    pub pwr_out_max: si::Power,
    /// starting/baseline transient power limit
    #[serde(default)]
    pub pwr_out_max_init: si::Power,
    // TODO: consider a ramp down rate, which may be needed for fuel cells
    /// lag time for ramp up
    pub pwr_ramp_lag: si::Time,
    /// interpolator for calculating [Self] efficiency as a function of output power
    #[api(skip_get, skip_set)]
    pub eff_interp_from_pwr_out: Interpolator,
    /// power at which peak efficiency occurs
    #[serde(skip)]
    pub pwr_for_peak_eff: si::Power,
    /// idle fuel power to overcome internal friction (not including aux load) \[W\]
    pub pwr_idle_fuel: si::Power,
    /// time step interval between saves. 1 is a good option. If None, no saving occurs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub save_interval: Option<usize>,
    /// struct for tracking current state
    #[serde(default)]
    #[serde(skip_serializing_if = "EqDefault::eq_default")]
    pub state: FuelConverterState,
    /// Custom vector of [Self::state]
    #[serde(default)]
    #[serde(skip_serializing_if = "FuelConverterStateHistoryVec::is_empty")]
    pub history: FuelConverterStateHistoryVec,
    #[serde(skip)]
    // phantom private field to prevent direct instantiation in other modules
    #[api(skip_get, skip_set)]
    pub(in super::super) _phantom: PhantomData<()>,
}

impl SetCumulative for FuelConverter {
    fn set_cumulative(&mut self, dt: si::Time) {
        self.state.set_cumulative(dt);
    }
}

impl SerdeAPI for FuelConverter {}
impl Init for FuelConverter {
    fn init(&mut self) -> anyhow::Result<()> {
        let _ = self.mass().with_context(|| anyhow!(format_dbg!()))?;
        self.thrml.init()?;
        self.state.init().with_context(|| anyhow!(format_dbg!()))?;
        let eff_max = self.eff_max()?;
        self.pwr_for_peak_eff = *self
            .eff_interp_from_pwr_out
            .x()
            .with_context(|| format_dbg!())?
            .get(
                self.eff_interp_from_pwr_out
                    .f_x()
                    .unwrap()
                    .iter()
                    .position(|&eff| eff * uc::R == eff_max)
                    .with_context(|| format_dbg!())?,
            )
            .with_context(|| format_dbg!())?
            * self.pwr_out_max;
        Ok(())
    }
}

impl Mass for FuelConverter {
    fn mass(&self) -> anyhow::Result<Option<si::Mass>> {
        let derived_mass = self
            .derived_mass()
            .with_context(|| anyhow!(format_dbg!()))?;
        if let (Some(derived_mass), Some(set_mass)) = (derived_mass, self.mass) {
            ensure!(
                utils::almost_eq_uom(&set_mass, &derived_mass, None),
                format!(
                    "{}",
                    format_dbg!(utils::almost_eq_uom(&set_mass, &derived_mass, None)),
                )
            );
        }
        Ok(self.mass)
    }

    fn set_mass(
        &mut self,
        new_mass: Option<si::Mass>,
        side_effect: MassSideEffect,
    ) -> anyhow::Result<()> {
        let derived_mass = self
            .derived_mass()
            .with_context(|| anyhow!(format_dbg!()))?;
        if let (Some(derived_mass), Some(new_mass)) = (derived_mass, new_mass) {
            if derived_mass != new_mass {
                match side_effect {
                    MassSideEffect::Extensive => {
                        self.pwr_out_max = self.specific_pwr.ok_or_else(|| {
                            anyhow!(
                                "{}\nExpected `self.specific_pwr` to be `Some`.",
                                format_dbg!()
                            )
                        })? * new_mass;
                    }
                    MassSideEffect::Intensive => {
                        self.specific_pwr = Some(self.pwr_out_max / new_mass);
                    }
                    MassSideEffect::None => {
                        self.specific_pwr = None;
                    }
                }
            }
        } else if new_mass.is_none() {
            self.specific_pwr = None;
        }
        self.mass = new_mass;
        Ok(())
    }

    fn derived_mass(&self) -> anyhow::Result<Option<si::Mass>> {
        Ok(self
            .specific_pwr
            .map(|specific_pwr| self.pwr_out_max / specific_pwr))
    }

    fn expunge_mass_fields(&mut self) {
        self.mass = None;
        self.specific_pwr = None;
    }
}

impl SaveInterval for FuelConverter {
    fn save_interval(&self) -> anyhow::Result<Option<usize>> {
        Ok(self.save_interval)
    }
    fn set_save_interval(&mut self, save_interval: Option<usize>) -> anyhow::Result<()> {
        self.save_interval = save_interval;
        Ok(())
    }
}

// non-py methods

impl FuelConverter {
    /// Sets maximum possible total power [FuelConverter]
    /// can produce.
    /// # Arguments
    /// - `dt`: time step size
    pub fn set_curr_pwr_out_max(&mut self, dt: si::Time) -> anyhow::Result<()> {
        if self.pwr_out_max_init == si::Power::ZERO {
            // TODO: think about how to initialize power
            self.pwr_out_max_init = self.pwr_out_max / 10.
        };
        self.state.pwr_out_max =
            (self.state.pwr_prop + self.state.pwr_aux + self.pwr_out_max / self.pwr_ramp_lag * dt)
                .min(self.pwr_out_max)
                .max(self.pwr_out_max_init);
        Ok(())
    }

    /// Sets maximum possible propulsion-related power [FuelConverter]
    /// can produce, accounting for any aux-related power required.
    /// # Arguments
    /// - `pwr_aux`: aux-related power required from this component
    pub fn set_curr_pwr_prop_max(&mut self, pwr_aux: si::Power) -> anyhow::Result<()> {
        ensure!(
            pwr_aux >= si::Power::ZERO,
            format!(
                "{}\n`pwr_aux` must be >= 0",
                format_dbg!(pwr_aux >= si::Power::ZERO),
            )
        );
        self.state.pwr_aux = pwr_aux;
        self.state.pwr_prop_max = self.state.pwr_out_max - pwr_aux;
        Ok(())
    }

    /// Solves for this powertrain system/component efficiency and sets/returns power output values.
    /// # Arguments
    /// - `pwr_out_req`: tractive power output required to achieve presribed speed
    /// - `fc_on`: whether component is actively running
    /// - `dt`: time step size
    pub fn solve(
        &mut self,
        pwr_out_req: si::Power,
        fc_on: bool,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        self.state.fc_on = fc_on;
        if fc_on {
            self.state.time_on += dt;
        } else {
            self.state.time_on = si::Time::ZERO;
        }
        // NOTE: think about the possibility of engine braking, not urgent
        ensure!(
            pwr_out_req >= si::Power::ZERO,
            format!(
                "{}\n`pwr_out_req` must be >= 0",
                format_dbg!(pwr_out_req >= si::Power::ZERO),
            )
        );
        // if the engine is not on, `pwr_out_req` should be 0.0
        ensure!(
            fc_on || (pwr_out_req == si::Power::ZERO && self.state.pwr_aux == si::Power::ZERO),
            format!(
                "{}\nEngine is off but pwr_out_req + pwr_aux is non-zero\n`pwr_out_req`: {} kW\n`self.state.pwr_aux`: {} kW",
                format_dbg!(
                    fc_on
                        || (pwr_out_req == si::Power::ZERO
                            && self.state.pwr_aux == si::Power::ZERO)
                ),
               pwr_out_req.get::<si::kilowatt>(),
               self.state.pwr_aux.get::<si::kilowatt>()
            )
        );
        self.state.pwr_prop = pwr_out_req;
        self.state.eff = if fc_on {
            uc::R
                * self
                    .eff_interp_from_pwr_out
                    .interpolate(&[
                        ((pwr_out_req + self.state.pwr_aux) / self.pwr_out_max).get::<si::ratio>()
                    ])
                    .with_context(|| {
                        anyhow!(
                            "{}\n failed to calculate {}",
                            format_dbg!(),
                            stringify!(self.state.eff)
                        )
                    })?
        } else {
            si::Ratio::ZERO
        };
        ensure!(
            (self.state.eff >= 0.0 * uc::R && self.state.eff <= 1.0 * uc::R),
            format!(
                "fc efficiency ({}) must be either between 0 and 1",
                self.state.eff.get::<si::ratio>()
            )
        );

        // TODO: consider how idle is handled.  The goal is to make it so that even if `self.state.pwr_aux` is
        // zero, there will be fuel consumption to overcome internal dissipation.
        self.state.pwr_fuel = if self.state.fc_on {
            ((pwr_out_req + self.state.pwr_aux) / self.state.eff).max(self.pwr_idle_fuel)
        } else {
            si::Power::ZERO
        };
        self.state.pwr_loss = self.state.pwr_fuel - self.state.pwr_prop;

        // TODO: put this in `SetCumulative::set_custom_cumulative`
        // ensure!(
        //     self.state.energy_loss.get::<si::joule>() >= 0.0,
        //     format!(
        //         "{}\nEnergy loss must be non-negative",
        //         format_dbg!(self.state.energy_loss.get::<si::joule>() >= 0.0)
        //     )
        // );
        Ok(())
    }

    pub fn solve_thermal(
        &mut self,
        te_amb: si::Temperature,
        heat_demand: si::Power,
        veh_state: VehicleState,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        let veh_speed = veh_state.speed_ach;
        self.thrml
            .solve(&self.state, te_amb, heat_demand, veh_speed, dt)
            .with_context(|| format_dbg!())?;
        Ok(())
    }

    pub fn eff_max(&self) -> anyhow::Result<si::Ratio> {
        Ok(self
            .eff_interp_from_pwr_out
            .f_x()
            .with_context(|| format_dbg!())?
            .iter()
            .fold(f64::NEG_INFINITY, |acc, &curr| acc.max(curr))
            * uc::R)
    }

    /// If thermal model is appropriately configured, returns current lumped engine temperature
    pub fn temperature(&self) -> Option<si::Temperature> {
        match &self.thrml {
            FuelConverterThermalOption::FuelConverterThermal(fct) => Some(fct.state.temperature),
            FuelConverterThermalOption::None => None,
        }
    }
}

// impl FuelConverter {
//     impl_get_set_eff_max_min!();
//     impl_get_set_eff_range!();
// }

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
pub struct FuelConverterState {
    /// time step index
    pub i: usize,
    /// max total output power fc can produce at current time
    pub pwr_out_max: si::Power,
    /// max propulsion power fc can produce at current time
    pub pwr_prop_max: si::Power,
    /// efficiency evaluated at current demand
    pub eff: si::Ratio,
    /// instantaneous power going to drivetrain, not including aux
    pub pwr_prop: si::Power,
    /// integral of [Self::pwr_prop]
    pub energy_prop: si::Energy,
    /// power going to auxiliaries
    pub pwr_aux: si::Power,
    /// Integral of [Self::pwr_aux]
    pub energy_aux: si::Energy,
    /// instantaneous fuel power flow
    pub pwr_fuel: si::Power,
    /// Integral of [Self::pwr_fuel]
    pub energy_fuel: si::Energy,
    /// loss power, including idle
    pub pwr_loss: si::Power,
    /// Integral of [Self::pwr_loss]
    pub energy_loss: si::Energy,
    /// If true, engine is on, and if false, off (no idle)
    pub fc_on: bool,
    /// Time the engine has been on
    pub time_on: si::Time,
}

impl SerdeAPI for FuelConverterState {}
impl Init for FuelConverterState {}

/// Options for handling [FuelConverter] thermal model
#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, IsVariant)]
pub enum FuelConverterThermalOption {
    /// Basic thermal plant for [FuelConverter]
    FuelConverterThermal(Box<FuelConverterThermal>),
    /// no thermal plant for [FuelConverter]
    #[default]
    None,
}
impl Init for FuelConverterThermalOption {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::FuelConverterThermal(fct) => fct.init()?,
            Self::None => {}
        }
        Ok(())
    }
}
impl SerdeAPI for FuelConverterThermalOption {}
impl FuelConverterThermalOption {
    /// Solve change in temperature and other thermal effects
    /// # Arguments
    /// - `fc_state`: [FuelConverter] state
    /// - `te_amb`: ambient temperature
    /// - `heat_demand`: heat demand from HVAC system
    /// - `veh_speed`: current achieved speed
    fn solve(
        &mut self,
        fc_state: &FuelConverterState,
        te_amb: si::Temperature,
        heat_demand: si::Power,
        veh_speed: si::Velocity,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        match self {
            Self::FuelConverterThermal(fct) => fct
                .solve(fc_state, te_amb, heat_demand, veh_speed, dt)
                .with_context(|| format_dbg!())?,
            Self::None => {}
        }
        Ok(())
    }
}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// Struct for modeling Fuel Converter (e.g. engine, fuel cell.)
pub struct FuelConverterThermal {
    /// [FuelConverter] thermal capacitance
    pub heat_capacitance: si::HeatCapacity,
    /// parameter for engine characteristic length for heat transfer calcs
    pub length_for_convection: si::Length,
    /// parameter for heat transfer coeff from [FuelConverter] to ambient during vehicle stop
    pub htc_to_amb_stop: si::HeatTransferCoeff,

    /// Heat transfer coefficient between adiabatic flame temperature and [FuelConverterThermal] temperature
    pub htc_from_comb: si::ThermalConductance,
    /// Max coefficient for fraction of combustion heat that goes to [FuelConverter]
    /// (engine) thermal mass. Remainder goes to environment (e.g. via tailpipe).
    pub max_frac_from_comb: si::Ratio,
    /// parameter for temperature at which thermostat starts to open
    #[api(skip_get, skip_set)]
    pub tstat_te_sto: Option<si::Temperature>,
    /// temperature delta over which thermostat is partially open
    #[api(skip_get, skip_set)]
    pub tstat_te_delta: Option<si::Temperature>,
    #[serde(skip_serializing, deserialize_with = "tstat_interp_default")]
    #[api(skip_get, skip_set)]
    pub tstat_interp: Interp1D,
    /// Radiator effectiveness -- ratio of active heat rejection from
    /// radiator to passive heat rejection, always greater than 1
    pub radiator_effectiveness: si::Ratio,
    /// struct for tracking current state
    #[serde(default)]
    #[serde(skip_serializing_if = "EqDefault::eq_default")]
    pub state: FuelConverterThermalState,
    /// Custom vector of [Self::state]
    #[serde(default)]
    #[serde(skip_serializing_if = "FuelConverterThermalStateHistoryVec::is_empty")]
    pub history: FuelConverterThermalStateHistoryVec,
}

/// Dummy interpolator that will be overridden in [FuelConverterThermal::init]
fn tstat_interp_default<'de, D>(_deserializer: D) -> Result<Interp1D, D::Error>
where
    D: Deserializer<'de>,
{
    Ok(Interp1D::new(
        vec![85.0, 90.0],
        vec![0.0, 1.0],
        Strategy::Linear,
        Extrapolate::Clamp,
    )
    .with_context(|| format_dbg!())
    .unwrap())
}

lazy_static! {
    /// gasoline stoichiometric air-fuel ratio https://en.wikipedia.org/wiki/Air%E2%80%93fuel_ratio
    pub static ref AFR_STOICH_GASOLINE: si::Ratio = uc::R * 14.7;
    /// gasoline density in https://inchem.org/documents/icsc/icsc/eics1400.htm
    /// This is reasonably constant with respect to temperature and pressure
    pub static ref GASOLINE_DENSITY: si::MassDensity = 0.75 * uc::KG / uc::L;
    /// TODO: find a source for this value
    pub static ref GASOLINE_LHV: si::SpecificEnergy = 33.7 * uc::KWH / uc::GALLON / *GASOLINE_DENSITY;
    pub static ref TE_ADIABATIC_STD: si::Temperature= Air::get_te_from_u(
            Air::get_specific_energy(*TE_STD_AIR).with_context(|| format_dbg!()).unwrap()
                + (Octane::get_specific_energy(*TE_STD_AIR).with_context(|| format_dbg!()).unwrap()
                    + *GASOLINE_LHV)
                    / *AFR_STOICH_GASOLINE,
        )
        .with_context(|| format_dbg!()).unwrap();
}

impl FuelConverterThermal {
    /// Solve change in temperature and other thermal effects
    /// # Arguments
    /// - `fc_state`: [FuelConverter] state
    /// - `te_amb`: ambient temperature
    /// - `heat_demand`: heat demand from HVAC system
    /// - `veh_speed`: current achieved speed
    /// - `dt`: simulation time step size
    fn solve(
        &mut self,
        fc_state: &FuelConverterState,
        te_amb: si::Temperature,
        heat_demand: si::Power,
        veh_speed: si::Velocity,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        // film temperature for external convection calculations
        let te_air_film = 0.5 * (self.state.temperature + te_amb);
        // Reynolds number = density * speed * diameter / dynamic viscosity
        // NOTE: might be good to pipe in elevation
        let fc_air_film_re =
            Air::get_density(Some(te_air_film), None) * veh_speed * self.length_for_convection
                / Air::get_dyn_visc(te_air_film).with_context(|| format_dbg!())?;

        // calculate heat transfer coeff. from engine to ambient [W / (m ** 2 * K)]
        self.state.htc_to_amb = if veh_speed < 1.0 * uc::MPS {
            // if stopped, scale based on thermostat opening and constant convection
            (uc::R
                + self
                    .tstat_interp
                    .interpolate(&[self.state.temperature.get::<si::degree_celsius>()])
                    .with_context(|| format_dbg!())?
                    * self.radiator_effectiveness)
                * self.htc_to_amb_stop
        } else {
            // Calculate heat transfer coefficient for sphere,
            // from Incropera's Intro to Heat Transfer, 5th Ed., eq. 7.44
            let sphere_conv_params = get_sphere_conv_params(fc_air_film_re.get::<si::ratio>());
            let htc_to_amb_sphere: si::HeatTransferCoeff = sphere_conv_params.0
                * fc_air_film_re.get::<si::ratio>().powf(sphere_conv_params.1)
                * Air::get_pr(te_air_film)
                    .with_context(|| format_dbg!())?
                    .get::<si::ratio>()
                    .powf(1.0 / 3.0)
                * Air::get_therm_cond(te_air_film).with_context(|| format_dbg!())?
                / self.length_for_convection;
            // if stopped, scale based on thermostat opening and constant convection
            self.tstat_interp
                .interpolate(&[self.state.temperature.get::<si::degree_celsius>()])
                .with_context(|| format_dbg!())?
                * htc_to_amb_sphere
        };

        self.state.heat_to_amb =
            self.state.htc_to_amb * PI * self.length_for_convection.powi(typenum::P2::new()) / 4.0
                * (self.state.temperature - te_amb);

        // let heat_to_amb = ;
        // assumes fuel/air mixture is entering combustion chamber at block temperature
        // assumes stoichiometric combustion
        self.state.te_adiabatic = Air::get_te_from_u(
            Air::get_specific_energy(self.state.temperature).with_context(|| format_dbg!())?
                + (Octane::get_specific_energy(self.state.temperature)
                    .with_context(|| format_dbg!())?
                    + *GASOLINE_LHV)
                    / *AFR_STOICH_GASOLINE,
        )
        .with_context(|| format_dbg!())?;
        // heat that will go both to the block and out the exhaust port
        let heat_gen = fc_state.pwr_fuel - fc_state.pwr_prop;
        let delta_temp: si::Temperature = (((self.htc_from_comb
            * (self.state.te_adiabatic - self.state.temperature))
            .min(self.max_frac_from_comb * heat_gen)
            - heat_demand
            - self.state.heat_to_amb)
            * dt)
            / self.heat_capacitance;
        self.state.temperature += delta_temp;
        Ok(())
    }
}
impl SerdeAPI for FuelConverterThermal {}
impl Init for FuelConverterThermal {
    fn init(&mut self) -> anyhow::Result<()> {
        self.tstat_te_sto = self
            .tstat_te_sto
            .or(Some(85. * uc::KELVIN + *uc::CELSIUS_TO_KELVIN));
        self.tstat_te_delta = self.tstat_te_delta.or(Some(5. * uc::KELVIN));
        self.tstat_interp = Interp1D::new(
            vec![
                self.tstat_te_sto.unwrap().get::<si::degree_celsius>(),
                self.tstat_te_sto.unwrap().get::<si::degree_celsius>()
                    + self.tstat_te_delta.unwrap().get::<si::degree_celsius>(),
            ],
            vec![0.0, 1.0],
            Strategy::Linear,
            Extrapolate::Clamp,
        )
        .with_context(|| format_dbg!())?;
        Ok(())
    }
}

#[fastsim_api]
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative)]
pub struct FuelConverterThermalState {
    /// time step index
    pub i: usize,
    /// Adiabatic flame temperature assuming complete (i.e. all fuel is consumed
    /// if fuel lean or stoich or all air is consumed if fuel rich) combustion
    pub te_adiabatic: si::Temperature,
    /// Current engine thermal mass temperature (lumped engine block and coolant)
    pub temperature: si::Temperature,
    /// Current heat transfer coefficient from [FuelConverter] to ambient
    pub htc_to_amb: si::HeatTransferCoeff,
    /// Current heat transfer to ambient
    pub heat_to_amb: si::Power,
}

impl Init for FuelConverterThermalState {}
impl SerdeAPI for FuelConverterThermalState {}
impl Default for FuelConverterThermalState {
    fn default() -> Self {
        Self {
            i: Default::default(),
            te_adiabatic: *TE_ADIABATIC_STD,
            temperature: *TE_STD_AIR,
            htc_to_amb: Default::default(),
            heat_to_amb: Default::default(),
        }
    }
}
