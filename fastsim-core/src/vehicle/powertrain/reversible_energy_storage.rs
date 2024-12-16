use super::*;

#[allow(unused_imports)]
#[cfg(feature = "pyo3")]
use crate::pyo3::*;

const TOL: f64 = 1e-3;

#[fastsim_api(
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

    #[getter("eff_range")]
    fn get_eff_range_py(&self) -> f64 {
        self.get_eff_range()
    }

    // #[setter("__eff_range")]
    // fn set_eff_range_py(&mut self, eff_range: f64) -> anyhow::Result<()> {
    //     self.set_eff_range(eff_range)
    // }

    // TODO: decide on way to deal with `side_effect` coming after optional arg and uncomment
    #[pyo3(name = "set_mass")]
    fn set_mass_py(&mut self, mass_kg: Option<f64>, side_effect: Option<String>) -> anyhow::Result<()> {
        let side_effect = side_effect.unwrap_or_else(|| "Intensive".into());
        self.set_mass(
            mass_kg.map(|m| m * uc::KG),
            MassSideEffect::try_from(side_effect)?
        )?;
        Ok(())
    }

    #[getter("mass_kg")]
    fn get_mass_kg_py(&mut self) -> anyhow::Result<Option<f64>> {
        Ok(self.mass()?.map(|m| m.get::<si::kilogram>()))
    }

    #[getter]
    fn get_specific_energy_kjoules_per_kg(&self) -> Option<f64> {
        self.specific_energy.map(|se| se.get::<si::kilojoule_per_kilogram>())
    }

    #[getter]
    fn get_energy_capacity_usable_joules(&self) -> f64 {
        self.energy_capacity_usable().get::<si::joule>()
    }

    #[pyo3(name = "set_default_1d_interp")]
    fn set_default_1d_interp_py(&mut self) -> anyhow::Result<()> {
        self.set_default_1d_interp()
    }

    #[pyo3(name = "set_default_2d_interp")]
    fn set_default_2d_interp_py(&mut self) -> anyhow::Result<()> {
        self.set_default_2d_interp()
    }

    #[pyo3(name = "set_default_3d_interp")]
    fn set_default_3d_interp_py(&mut self) -> anyhow::Result<()> {
        self.set_default_3d_interp()
    }
)]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
#[non_exhaustive]
/// Struct for modeling technology-naive Reversible Energy Storage (e.g. battery, flywheel).
pub struct ReversibleEnergyStorage {
    /// [Self] Thermal plant, including thermal management controls
    #[serde(default, skip_serializing_if = "RESThermalOption::is_none")]
    #[api(skip_get, skip_set)]
    pub thrml: RESThermalOption,
    /// ReversibleEnergyStorage mass
    #[serde(default)]
    #[api(skip_get, skip_set)]
    pub(in super::super) mass: Option<si::Mass>,
    /// ReversibleEnergyStorage specific energy
    #[api(skip_get, skip_set)]
    pub(in super::super) specific_energy: Option<si::SpecificEnergy>,
    /// Max output (and input) power battery can produce (accept)
    pub pwr_out_max: si::Power,

    /// Total energy capacity of battery of full discharge SOC of 0.0 and 1.0
    pub energy_capacity: si::Energy,

    /// interpolator for calculating [Self] efficiency as a function of the following variants:  
    /// - 0d -- constant -- handled on a round trip basis
    /// - 1d -- linear w.r.t. power
    /// - 2d -- linear w.r.t. power and SOC
    /// - 3d -- linear w.r.t. power, SOC, and temperature
    #[api(skip_get, skip_set)]
    pub eff_interp: Interpolator,

    /// Hard limit on minimum SOC, e.g. 0.05
    pub min_soc: si::Ratio,
    /// Hard limit on maximum SOC, e.g. 0.95
    pub max_soc: si::Ratio,

    /// Time step interval at which history is saved
    #[serde(skip_serializing_if = "Option::is_none")]
    pub save_interval: Option<usize>,
    /// struct for tracking current state
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: ReversibleEnergyStorageState,
    /// Custom vector of [Self::state]
    #[serde(
        default,
        skip_serializing_if = "ReversibleEnergyStorageStateHistoryVec::is_empty"
    )]
    pub history: ReversibleEnergyStorageStateHistoryVec,
}

impl ReversibleEnergyStorage {
    pub fn solve(&mut self, pwr_out_req: si::Power, dt: si::Time) -> anyhow::Result<()> {
        let te_res = self.temperature();
        let state = &mut self.state;

        ensure!(
            state.soc <= self.max_soc || (pwr_out_req + state.pwr_aux) >= si::Power::ZERO,
            "{}\n{}",
            format_dbg!(pwr_out_req + state.pwr_aux),
            state.soc.get::<si::ratio>()
        );
        ensure!(
            state.soc >= self.min_soc || (pwr_out_req + state.pwr_aux) <= si::Power::ZERO,
            "{}\n{}",
            format_dbg!(pwr_out_req + state.pwr_aux),
            state.soc.get::<si::ratio>()
        );

        state.pwr_out_prop = pwr_out_req;
        state.pwr_out_electrical = state.pwr_out_prop + state.pwr_aux;

        if pwr_out_req + state.pwr_aux >= si::Power::ZERO {
            // discharging
            ensure!(
                utils::almost_le_uom(&(pwr_out_req + state.pwr_aux), &self.pwr_out_max, Some(TOL)),
                "{}\nres required power ({:.6} kW) exceeds static max discharge power ({:.6} kW)\nstate.soc = {}",
                format_dbg!(utils::almost_le_uom(
                    &(pwr_out_req + state.pwr_aux),
                    &self.pwr_out_max,
                    Some(TOL)
                )),
                (pwr_out_req + state.pwr_aux).get::<si::kilowatt>(),
                state.pwr_disch_max.get::<si::kilowatt>(),
                state.soc.get::<si::ratio>()
            );
            ensure!(
                utils::almost_le_uom(&(pwr_out_req + state.pwr_aux), &state.pwr_disch_max, Some(TOL)),
                "{}\nres required power ({:.6} kW) exceeds current max discharge power ({:.6} kW)\nstate.soc = {}",
                format_dbg!(utils::almost_le_uom(&(pwr_out_req + state.pwr_aux), &state.pwr_disch_max, Some(TOL))),
                (pwr_out_req + state.pwr_aux).get::<si::kilowatt>(),
                state.pwr_disch_max.get::<si::kilowatt>(),
                state.soc.get::<si::ratio>()
            );
        } else {
            // charging
            ensure!(
                utils::almost_ge_uom(
                    &(pwr_out_req + state.pwr_aux),
                    &-self.pwr_out_max,
                    Some(TOL)
                ),
                format!(
                    "{}\nres required power ({:.6} kW) exceeds static max power ({:.6} kW)",
                    format_dbg!(utils::almost_ge_uom(
                        &(pwr_out_req + state.pwr_aux),
                        &-self.pwr_out_max,
                        Some(TOL)
                    )),
                    (pwr_out_req + state.pwr_aux).get::<si::kilowatt>(),
                    state.pwr_charge_max.get::<si::kilowatt>()
                )
            );
            ensure!(
                utils::almost_ge_uom(
                    &(pwr_out_req + state.pwr_aux),
                    &-state.pwr_charge_max,
                    Some(TOL)
                ),
                format!(
                    "{}\nres required power ({:.6} kW) exceeds current max power ({:.6} kW)",
                    format_dbg!(utils::almost_ge_uom(
                        &(pwr_out_req + state.pwr_aux),
                        &-state.pwr_charge_max,
                        Some(TOL)
                    )),
                    (pwr_out_req + state.pwr_aux).get::<si::kilowatt>(),
                    state.pwr_charge_max.get::<si::kilowatt>()
                )
            );
        }
        state.eff = match &self.eff_interp {
            Interpolator::Interp0D(eff) => *eff * uc::R,
            Interpolator::Interp1D(interp1d) => {
                interp1d.interpolate(&[state.pwr_out_electrical.get::<si::watt>()])? * uc::R
            }
            Interpolator::Interp2D(interp2d) => {
                interp2d.interpolate(&[
                    state.pwr_out_electrical.get::<si::watt>(),
                    state.soc.get::<si::ratio>(),
                ])? * uc::R
            }
            Interpolator::Interp3D(interp3d) => {
                interp3d.interpolate(&[
                    state.pwr_out_electrical.get::<si::watt>(),
                    state.soc.get::<si::ratio>(),
                    te_res
                        .with_context(|| format_dbg!("Expected thermal model to be configured"))?
                        .get::<si::degree_celsius>(),
                ])? * uc::R
            }
            _ => bail!("Invalid interpolator.  See docs for `ReversibleEnergyStorage::eff_interp`"),
        };
        ensure!(
            state.eff >= 0.0 * uc::R && state.eff <= 1.0 * uc::R,
            format!(
                "{}\nres efficiency ({}) must be between 0 and 1",
                format_dbg!(state.eff >= 0.0 * uc::R && state.eff <= 1.0 * uc::R),
                state.eff.get::<si::ratio>()
            )
        );

        state.pwr_out_chemical = if state.pwr_out_electrical > si::Power::ZERO {
            // if positive, chemical power must be greater than electrical power
            // i.e. not all chemical power can be converted to electrical power
            state.pwr_out_electrical / state.eff
        } else {
            // if negative, chemical power, must be less than electrical power
            // i.e. not all electrical power can be converted back to chemical power
            state.pwr_out_electrical * state.eff
        };

        state.pwr_loss = (state.pwr_out_chemical - state.pwr_out_electrical).abs();

        state.soc -= state.pwr_out_chemical * dt / self.energy_capacity;

        Ok(())
    }

    /// Solve change in temperature and other thermal effects
    /// # Arguments
    /// - `fc_state`: [ReversibleEnergyStorage] state
    /// - `te_amb`: ambient temperature
    /// - `dt`: time step size
    pub fn solve_thermal(&mut self, te_amb: si::Temperature, dt: si::Time) -> anyhow::Result<()> {
        self.thrml
            .solve(self.state, te_amb, dt)
            .with_context(|| format_dbg!())
    }

    /// Sets and returns max output and max regen power based on current state
    /// # Arguments
    /// - `dt`: time step size
    /// - `disch_buffer`: buffer offset from static SOC limit at which discharging is not allowed
    /// - `chrg_buffer`: buffer offset from static SOC limit at which charging is not allowed
    pub fn set_curr_pwr_out_max(
        &mut self,
        dt: si::Time,
        disch_buffer: si::Energy,
        chrg_buffer: si::Energy,
    ) -> anyhow::Result<()> {
        self.set_pwr_disch_max(dt, disch_buffer)?;
        self.set_pwr_charge_max(dt, chrg_buffer)?;

        Ok(())
    }

    /// # Arguments
    /// - `dt`: time step size
    /// - `buffer`: buffer below static maximum SOC above which charging is disabled
    pub fn set_pwr_charge_max(
        &mut self,
        dt: si::Time,
        chrg_buffer: si::Energy,
    ) -> anyhow::Result<()> {
        // to protect against excessive topping off of the battery
        let soc_buffer_delta = (chrg_buffer
            / (self.energy_capacity * (self.max_soc - self.min_soc)))
            .max(si::Ratio::ZERO);
        ensure!(soc_buffer_delta >= si::Ratio::ZERO, "{}", format_dbg!());
        self.state.soc_regen_buffer = self.max_soc - soc_buffer_delta;
        let pwr_max_for_dt =
            ((self.max_soc - self.state.soc) * self.energy_capacity / dt).max(si::Power::ZERO);
        self.state.pwr_charge_max = if self.state.soc <= self.state.soc_regen_buffer {
            self.pwr_out_max
        } else if self.state.soc < self.max_soc && soc_buffer_delta > si::Ratio::ZERO {
            self.pwr_out_max * (self.max_soc - self.state.soc) / soc_buffer_delta
        } else {
            // current SOC is less than both
            si::Power::ZERO
        }
        .min(pwr_max_for_dt);

        ensure!(
            self.state.pwr_charge_max >= si::Power::ZERO,
            "{}\n`{}` ({} W) must be greater than or equal to zero\n{}",
            format_dbg!(),
            stringify!(self.state.pwr_charge_max),
            self.state.pwr_charge_max.get::<si::watt>().format_eng(None),
            format_dbg!(soc_buffer_delta)
        );

        Ok(())
    }

    /// # Arguments
    /// - `dt`: time step size
    /// - `buffer`: buffer above static minimum SOC above which charging is disabled
    pub fn set_pwr_disch_max(
        &mut self,
        dt: si::Time,
        disch_buffer: si::Energy,
    ) -> anyhow::Result<()> {
        // to protect against excessive bottoming out of the battery
        let soc_buffer_delta = (disch_buffer / self.energy_capacity_usable()).max(si::Ratio::ZERO);
        ensure!(soc_buffer_delta >= si::Ratio::ZERO, "{}", format_dbg!());
        self.state.soc_disch_buffer = self.min_soc + soc_buffer_delta;
        let pwr_max_for_dt =
            ((self.state.soc - self.min_soc) * self.energy_capacity / dt).max(si::Power::ZERO);
        self.state.pwr_disch_max = if self.state.soc > self.state.soc_disch_buffer {
            self.pwr_out_max
        } else if self.state.soc > self.min_soc && soc_buffer_delta > si::Ratio::ZERO {
            self.pwr_out_max * (self.state.soc - self.min_soc) / soc_buffer_delta
        } else {
            // current SOC is less than both
            si::Power::ZERO
        }
        .min(pwr_max_for_dt);

        ensure!(
            self.state.pwr_disch_max >= si::Power::ZERO,
            "{}\n`{}` ({} W) must be greater than or equal to zero\n{}",
            format_dbg!(),
            stringify!(self.state.pwr_disch_max),
            self.state.pwr_disch_max.get::<si::watt>().format_eng(None),
            format_dbg!(soc_buffer_delta)
        );

        Ok(())
    }

    /// Set current maximum power available for propulsion
    /// # Arguments
    /// - `pwr_aux`: aux power demand on `ReversibleEnergyStorage`
    pub fn set_curr_pwr_prop_max(&mut self, pwr_aux: si::Power) -> anyhow::Result<()> {
        let state = &mut self.state;
        state.pwr_aux = pwr_aux;
        state.pwr_prop_max = state.pwr_disch_max - pwr_aux;
        state.pwr_regen_max = state.pwr_charge_max + pwr_aux;

        ensure!(
            pwr_aux <= state.pwr_disch_max,
            "`{}` ({} W) must always be less than or equal to {} ({} W)\nsoc:{}",
            stringify!(pwr_aux),
            pwr_aux.get::<si::watt>().format_eng(None),
            stringify!(state.pwr_disch_max),
            state.pwr_disch_max.get::<si::watt>().format_eng(None),
            state.soc.get::<si::ratio>()
        );
        ensure!(
            state.pwr_prop_max >= si::Power::ZERO,
            "`{}` ({} W) must be greater than or equal to zero",
            stringify!(state.pwr_prop_max),
            state.pwr_prop_max.get::<si::watt>().format_eng(None)
        );
        ensure!(
            state.pwr_regen_max >= si::Power::ZERO,
            "`{}` ({} W) must be greater than or equal to zero",
            stringify!(state.pwr_regen_max),
            state.pwr_regen_max.get::<si::watt>().format_eng(None)
        );

        Ok(())
    }

    /// Sets specific energy and either mass or energy capacity of battery
    /// # Arguments
    /// - `specific_energy`: specific energy of battery
    /// - `side_effect`: whether to update mass or energy capacity
    pub fn set_specific_energy(
        mut self,
        specific_energy: si::SpecificEnergy,
        side_effect: SpecificEnergySideEffect,
    ) -> anyhow::Result<()> {
        self.specific_energy = Some(specific_energy);
        match side_effect {
            SpecificEnergySideEffect::Mass => self.set_mass(
                Some(self.energy_capacity / specific_energy),
                MassSideEffect::Intensive,
            )?,
            SpecificEnergySideEffect::Energy => {
                self.energy_capacity = specific_energy
                    * self.mass.with_context(|| {
                        format_dbg!("Expected `ReversibleEnergyStorage::mass` to have been set.")
                    })?;
            }
        }
        Ok(())
    }

    pub fn get_eff_max(&self) -> f64 {
        todo!("adapt from ALTRIOS");
    }

    /// Scales eff_interp by ratio of new `eff_max` per current calculated
    /// max linearly, such that `eff_min` is untouched
    pub fn set_eff_max(&mut self, _eff_max: f64) -> Result<(), String> {
        todo!("adapt from ALTRIOS");
    }

    pub fn get_eff_min(&self) -> f64 {
        todo!("adapt from ALTRIOS");
    }

    /// Max value of `eff_interp` minus min value of `eff_interp`.
    pub fn get_eff_range(&self) -> f64 {
        self.get_eff_max() - self.get_eff_min()
    }

    /// Scales values of `eff_interp` without changing max such that max - min
    /// is equal to new range
    pub fn set_eff_range(&mut self, _eff_range: f64) -> anyhow::Result<()> {
        todo!("adapt from ALTRIOS");
    }

    /// Usable energy capacity, accounting for SOC limits
    pub fn energy_capacity_usable(&self) -> si::Energy {
        self.energy_capacity * (self.max_soc - self.min_soc)
    }

    /// Sets the ReversibleEnergyStorage eff_interp Interpolator to be a 1D
    /// interpolator with the default x and f_x arrays  
    /// Source of default efficiency values:  
    /// x: values in the third sub-array (corresponding to power) in Altrios'
    /// eta_interp_grid  
    /// f_x: efficiency array as a function of power at constant 50% SOC and 23
    /// degrees C corresponds to eta_interp_values[0][5] in Altrios
    #[cfg(all(feature = "yaml", feature = "resources"))]
    pub fn set_default_1d_interp(&mut self) -> anyhow::Result<()> {
        self.eff_interp = ninterp::Interpolator::from_resource("res/default_1d.yaml", false)?;
        Ok(())
    }

    /// Sets the ReversibleEnergyStorage eff_interp Interpolator to be a 2D
    /// interpolator with the default x, y and f_xy arrays  
    /// Source of default efficiency values:  
    /// x: values in the third sub-array (corresponding to power) in Altrios'
    /// eta_interp_grid  
    /// y: values in the second sub-array (corresponding to SOC) in
    /// Altrios' eta_interp_grid  
    /// f_xy: efficiency array as a function of power and SOC at constant 23
    /// degrees C corresponds to eta_interp_values[0] in Altrios, transposed so
    /// that the outermost layer is now power, and the innermost layer SOC (in
    /// altrios, the outermost layer is SOC and innermost is power)
    #[cfg(all(feature = "yaml", feature = "resources"))]
    pub fn set_default_2d_interp(&mut self) -> anyhow::Result<()> {
        self.eff_interp = ninterp::Interpolator::from_resource("res/default_2d.yaml", false)?;
        Ok(())
    }

    /// Sets the ReversibleEnergyStorage eff_interp Interpolator to be a 3D
    /// interpolator with the default x, y, z and f_xyz arrays  
    /// Source of default efficiency values:  
    /// x: values in the third sub-array (corresponding to power) in Altrios'
    /// eta_interp_grid  
    /// y: values in the second sub-array (corresponding to SOC) in Altrios'
    /// eta_interp_grid  
    /// z: values in the first sub-array (corresponding to temperature) in
    /// Altrios' eta_interp_grid  
    /// f_xyz: efficiency array as a function of power, SOC, and temperature
    /// corresponds to eta_interp_values in Altrios, transposed so that the
    /// outermost layer is now power, and the innermost layer temperature (in
    /// altrios, the outermost layer is temperature and innermost is power)
    #[cfg(all(feature = "yaml", feature = "resources"))]
    pub fn set_default_3d_interp(&mut self) -> anyhow::Result<()> {
        self.eff_interp = ninterp::Interpolator::from_resource("res/default_3d.yaml", false)?;
        Ok(())
    }

    /// If thermal model is appropriately configured, returns current lumped [Self] temperature
    pub fn temperature(&self) -> Option<si::Temperature> {
        match &self.thrml {
            RESThermalOption::RESLumpedThermal(rest) => Some(rest.state.temperature),
            RESThermalOption::None => None,
        }
    }
}

impl SetCumulative for ReversibleEnergyStorage {
    fn set_cumulative(&mut self, dt: si::Time) {
        self.state.set_cumulative(dt);
    }
}

impl Mass for ReversibleEnergyStorage {
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
                        self.energy_capacity = self.specific_energy.ok_or_else(|| {
                            anyhow!(
                                "{}\nExpected `self.specific_energy` to be `Some`.",
                                format_dbg!()
                            )
                        })? * new_mass;
                    }
                    MassSideEffect::Intensive => {
                        self.specific_energy = Some(self.energy_capacity / new_mass);
                    }
                    MassSideEffect::None => {
                        self.specific_energy = None;
                    }
                }
            }
        } else if new_mass.is_none() {
            self.specific_energy = None;
        }
        self.mass = new_mass;

        Ok(())
    }

    fn derived_mass(&self) -> anyhow::Result<Option<si::Mass>> {
        Ok(self
            .specific_energy
            .map(|specific_energy| self.energy_capacity / specific_energy))
    }

    fn expunge_mass_fields(&mut self) {
        self.mass = None;
        self.specific_energy = None;
    }
}

impl SerdeAPI for ReversibleEnergyStorage {}
impl Init for ReversibleEnergyStorage {
    fn init(&mut self) -> anyhow::Result<()> {
        let _ = self.mass().with_context(|| anyhow!(format_dbg!()))?;
        self.state.init().with_context(|| anyhow!(format_dbg!()))?;
        // TODO: make some kind of data validation framework to replace this code.
        ensure!(
            self.max_soc > self.min_soc,
            format!(
                "{}\n`max_soc`: {} must be greater than `min_soc`: {}`",
                format_dbg!(),
                self.max_soc.get::<si::ratio>(),
                self.min_soc.get::<si::ratio>(),
            )
        );
        Ok(())
    }
}

/// Controls which parameter to update when setting specific energy
pub enum SpecificEnergySideEffect {
    /// update mass
    Mass,
    /// update energy
    Energy,
}

#[fastsim_api]
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative)]
#[non_exhaustive]
// component limits
/// ReversibleEnergyStorage state variables
pub struct ReversibleEnergyStorageState {
    // limits
    /// max output power for propulsion during positive traction
    pub pwr_prop_max: si::Power,
    /// max regen power for propulsion during negative traction
    pub pwr_regen_max: si::Power,
    /// max discharge power total
    pub pwr_disch_max: si::Power,
    /// max charge power on the output side
    pub pwr_charge_max: si::Power,

    /// time step index
    pub i: usize,

    /// state of charge (SOC)
    pub soc: si::Ratio,
    /// SOC at which [ReversibleEnergyStorage] regen power begins linearly
    /// derating as it approaches maximum SOC
    pub soc_regen_buffer: si::Ratio,
    /// SOC at which [ReversibleEnergyStorage] discharge power begins linearly
    /// derating as it approaches minimum SOC
    pub soc_disch_buffer: si::Ratio,
    /// Chemical <-> Electrical conversion efficiency based on current power demand
    pub eff: si::Ratio,
    /// State of Health (SOH)
    pub soh: f64,

    // TODO: add `pwr_out_neg_electrical` and `pwr_out_pos_electrical` and corresponding energies
    // powers to separately pin negative- and positive-power operation
    /// total electrical power; positive is discharging
    pub pwr_out_electrical: si::Power,
    /// electrical power going to propulsion
    pub pwr_out_prop: si::Power,
    /// electrical power going to aux loads
    pub pwr_aux: si::Power,
    /// power dissipated as loss
    pub pwr_loss: si::Power,
    /// chemical power; positive is discharging
    pub pwr_out_chemical: si::Power,

    // cumulative energies
    /// cumulative total electrical energy; positive is discharging
    pub energy_out_electrical: si::Energy,
    /// cumulative electrical energy going to propulsion
    pub energy_out_prop: si::Energy,
    /// cumulative electrical energy going to aux loads
    pub energy_aux: si::Energy,
    /// cumulative energy dissipated as loss
    pub energy_loss: si::Energy,
    /// cumulative chemical energy; positive is discharging
    pub energy_out_chemical: si::Energy,
}

impl Default for ReversibleEnergyStorageState {
    fn default() -> Self {
        Self {
            pwr_prop_max: si::Power::ZERO,
            pwr_regen_max: si::Power::ZERO,
            pwr_disch_max: si::Power::ZERO,
            pwr_charge_max: si::Power::ZERO,
            i: Default::default(),
            soc: uc::R * 0.5,
            soc_regen_buffer: uc::R * 1.,
            soc_disch_buffer: si::Ratio::ZERO,
            eff: si::Ratio::ZERO,
            soh: 0.,
            pwr_out_electrical: si::Power::ZERO,
            pwr_out_prop: si::Power::ZERO,
            pwr_aux: si::Power::ZERO,
            pwr_loss: si::Power::ZERO,
            pwr_out_chemical: si::Power::ZERO,
            energy_out_electrical: si::Energy::ZERO,
            energy_out_prop: si::Energy::ZERO,
            energy_aux: si::Energy::ZERO,
            energy_loss: si::Energy::ZERO,
            energy_out_chemical: si::Energy::ZERO,
        }
    }
}

impl Init for ReversibleEnergyStorageState {}
impl SerdeAPI for ReversibleEnergyStorageState {}

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, IsVariant)]
pub enum RESThermalOption {
    /// Basic thermal plant for [ReversibleEnergyStorage]
    RESLumpedThermal(Box<RESLumpedThermal>),
    /// no thermal plant for [ReversibleEnergyStorage]
    #[default]
    None,
}
impl Init for RESThermalOption {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::RESLumpedThermal(rest) => rest.init()?,
            Self::None => {}
        }
        Ok(())
    }
}
impl SerdeAPI for RESThermalOption {}
impl RESThermalOption {
    /// Solve change in temperature and other thermal effects
    /// # Arguments
    /// - `res_state`: [ReversibleEnergyStorage] state
    /// - `te_amb`: ambient temperature
    /// - `dt`: time step size
    fn solve(
        &mut self,
        res_state: ReversibleEnergyStorageState,
        te_amb: si::Temperature,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        match self {
            Self::RESLumpedThermal(rest) => rest
                .solve(res_state, te_amb, dt)
                .with_context(|| format_dbg!())?,
            Self::None => {
                // TODO: make sure this triggers error if appropriate
            }
        }
        Ok(())
    }
}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// Struct for modeling [ReversibleEnergyStorage] (e.g. battery) thermal plant
pub struct RESLumpedThermal {
    /// [ReversibleEnergyStorage] thermal capacitance
    pub heat_capacitance: si::HeatCapacity,
    /// parameter for heat transfer coeff from [ReversibleEnergyStorage::thrml] to ambient
    pub htc_to_amb: si::HeatTransferCoeff,
    /// parameter for heat transfer coeff from [ReversibleEnergyStorage::thrml] to cabin
    pub htc_to_cab: si::HeatTransferCoeff,
    /// Thermal management system
    pub cntrl_sys: RESThermalControlSystem,
    /// current state
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: RESThermalState,
    /// history of state
    #[serde(default, skip_serializing_if = "RESThermalStateHistoryVec::is_empty")]
    pub history: RESThermalStateHistoryVec,
    // TODO: add `save_interval` and associated methods
}

impl SerdeAPI for RESLumpedThermal {}
impl Init for RESLumpedThermal {}
impl RESLumpedThermal {
    fn solve(
        &mut self,
        res_state: ReversibleEnergyStorageState,
        te_amb: si::Temperature,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        todo!();
        Ok(())
    }
}

#[fastsim_api]
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative)]
pub struct RESThermalState {
    /// time step index
    pub i: usize,
    /// Current thermal mass temperature
    pub temperature: si::Temperature,
    /// Thermal mass temperature at previous time step
    pub temp_prev: si::Temperature,
}

impl Init for RESThermalState {}
impl SerdeAPI for RESThermalState {}
impl Default for RESThermalState {
    fn default() -> Self {
        Self {
            i: Default::default(),
            temperature: *TE_STD_AIR,
            temp_prev: *TE_STD_AIR,
        }
    }
}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
/// HVAC system for [LumpedCabin]
pub struct RESThermalControlSystem {
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
    pub heat_source: RESHeatSource,
    /// max allowed aux load
    pub pwr_aux_max: si::Power,
    /// coefficient of performance of vapor compression cycle
    #[serde(default, skip_serializing_if = "EqDefault::eq_default")]
    pub state: RESThermalControlSystemState,
    #[serde(
        default,
        skip_serializing_if = "RESThermalControlSystemStateHistoryVec::is_empty"
    )]
    pub history: RESThermalControlSystemStateHistoryVec,
    // TODO: add `save_interval` and associated methods
}
impl Init for RESThermalControlSystem {}
impl SerdeAPI for RESThermalControlSystem {}
impl RESThermalControlSystem {
    pub fn get_pwr_thermal(
        &mut self,
        te_amb_air: si::Temperature,
        res_thrml_state: RESThermalState,
        dt: si::Time,
    ) -> anyhow::Result<si::Power> {
        let pwr_from_hvac = if res_thrml_state.temperature <= self.te_set + self.te_deadband
            && res_thrml_state.temperature >= self.te_set - self.te_deadband
        {
            // inside deadband; no hvac power is needed

            self.state.pwr_i = si::Power::ZERO; // reset to 0.0
            self.state.pwr_p = si::Power::ZERO;
            self.state.pwr_d = si::Power::ZERO;
            si::Power::ZERO
        } else {
            let te_delta_vs_set = res_thrml_state.temperature - self.te_set;
            let te_delta_vs_amb: si::Temperature = res_thrml_state.temperature - te_amb_air;

            self.state.pwr_p = -self.p * te_delta_vs_set;
            self.state.pwr_i -= self.i * uc::W / uc::KELVIN / uc::S * te_delta_vs_set * dt;
            self.state.pwr_i = self.state.pwr_i.max(-self.pwr_i_max).min(self.pwr_i_max);
            self.state.pwr_d = -self.d * uc::J / uc::KELVIN
                * ((res_thrml_state.temperature - res_thrml_state.temp_prev) / dt);

            // https://en.wikipedia.org/wiki/Coefficient_of_performance#Theoretical_performance_limits
            // cop_ideal is t_h / (t_h - t_c) for heating
            // cop_ideal is t_c / (t_h - t_c) for cooling

            // divide-by-zero protection and realistic limit on COP
            let cop_ideal = if te_delta_vs_amb.abs() < 5.0 * uc::KELVIN {
                // cabin is cooler than ambient + threshold
                // TODO: make this `5.0` not hardcoded
                res_thrml_state.temperature / (5.0 * uc::KELVIN)
            } else {
                res_thrml_state.temperature / te_delta_vs_amb.abs()
            };
            self.state.cop = cop_ideal * self.frac_of_ideal_cop;
            assert!(self.state.cop > 0.0 * uc::R);

            if res_thrml_state.temperature > self.te_set + self.te_deadband {
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
                match &self.heat_source {
                    RESHeatSource::SelfHeating => {
                        todo!()
                    }
                    RESHeatSource::HeatPump => {
                        todo!()
                    }
                    RESHeatSource::None => {}
                }
                pwr_thermal_from_hvac
            }
        };
        Ok(pwr_from_hvac)
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub enum RESHeatSource {
    /// Self heating
    SelfHeating,
    /// Heat pump provides heat for HVAC system
    HeatPump,
    /// No active heat source for [RESThermal]
    None,
}
impl Init for RESHeatSource {}
impl SerdeAPI for RESHeatSource {}

#[fastsim_api]
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, HistoryVec, SetCumulative,
)]
pub struct RESThermalControlSystemState {
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
    /// Cumulative energy demand by HVAC system from thermal component (e.g. [FuelConverter])
    pub energy_thermal_req: si::Energy,
}
impl Init for RESThermalControlSystemState {}
impl SerdeAPI for RESThermalControlSystemState {}
