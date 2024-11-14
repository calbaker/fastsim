use super::imports::*;
use super::*;

lazy_static! {
    /// room temperature
    pub static ref TE_STD_AIR: si::ThermodynamicTemperature = (22. + 273.15) * uc::KELVIN;
    /// pressure of air at 180 m and 22 C
    pub static ref STD_PRESSURE_AIR: si::Pressure = 99_346.3 * uc::PASCAL;
    /// density of air at 180 m ASL and 22 C
    pub static ref STD_DENSITY_AIR: si::MassDensity = 1.172 * uc::KGPM3;
    /// ideal gas constant for air
    pub static ref R_AIR: si::SpecificHeatCapacity = 287.0 * uc::J_PER_KG_K;
    /// standard elevation above sea level
    pub static ref H_STD: si::Length = 180.0 * uc::M;
}

#[fastsim_api(
    /// Returns density of air [kg/m^3]
    /// Source: <https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmet.html>  
    ///
    /// # Equations used
    /// T = 15.04 - .00649 * h  
    /// p = 101.29 * [(T + 273.1)/288.08]^5.256  
    ///
    /// # Arguments  
    /// * `te_air_deg_c` - optional ambient temperature [°C] of air, defaults to 22 C
    /// * `h_m` - optional elevation [m] above sea level, defaults to 180 m
    #[staticmethod]
    #[pyo3(name = "get_density")]
    pub fn get_density_py(te_air_deg_c: Option<f64>, h_m: Option<f64>) -> f64 {
        Self::get_density(
            te_air_deg_c.map(|te_air_deg_c| (te_air_deg_c + 273.15) * uc::KELVIN),
            h_m.map(|h_m| h_m * uc::M),
        )
        .get::<si::kilogram_per_cubic_meter>()
    }

    /// Returns thermal conductivity [W/(m*K)] of air
    /// # Arguments
    /// - `te_air`: temperature [°C] of air
    #[pyo3(name = "get_therm_cond")]
    #[staticmethod]
    pub fn get_therm_cond_py(te_air: f64) -> anyhow::Result<f64> {
        Ok(Self::get_therm_cond(te_air * uc::KELVIN - *uc::CELSIUS_TO_KELVIN)?.get::<si::watt_per_meter_kelvin>())
    }

    /// Returns constant pressure specific heat [J/(kg*K)] of air
    /// # Arguments
    /// - `te_air`: temperature [°C] of air
    #[pyo3(name = "get_specific_heat_cp")]
    #[staticmethod]
    pub fn get_specific_heat_cp_py(te_air: f64) -> anyhow::Result<f64> {
        Ok(Self::get_specific_heat_cp(te_air * uc::KELVIN - *uc::CELSIUS_TO_KELVIN)?.get::<si::joule_per_kilogram_kelvin>())
    }

    /// Returns specific enthalpy [J/kg] of air  
    /// # Arguments  
    /// - `te_air`: temperature [°C] of air
    #[pyo3(name = "get_specific_enthalpy")]
    #[staticmethod]
    pub fn get_specific_enthalpy_py(te_air: f64) -> anyhow::Result<f64> {
        Ok(Self::get_specific_enthalpy(te_air * uc::KELVIN - *uc::CELSIUS_TO_KELVIN)?.get::<si::joule_per_kilogram>())
    }

    /// Returns thermal Prandtl number of air
    /// # Arguments
    /// - `te_air`: temperature [°C] of air     
    #[pyo3(name = "get_pr")]
    #[staticmethod]
    pub fn get_pr_py(te_air: f64) -> anyhow::Result<f64> {
        Ok(Self::get_pr(te_air * uc::KELVIN - *uc::CELSIUS_TO_KELVIN)?.get::<si::ratio>())
    }

    /// Returns dynamic viscosity \[Pa*s\] of air
    /// # Arguments
    /// te_air: temperature [°C] of air
    #[pyo3(name = "get_dyn_visc")]
    #[staticmethod]
    pub fn get_dyn_visc_py(te_air: f64) -> anyhow::Result<f64> {
        Ok(Self::get_dyn_visc(te_air * uc::KELVIN - *uc::CELSIUS_TO_KELVIN)?.get::<si::pascal_second>())
    }

    /// Returns temperature [°C] of air
    /// # Arguments
    /// - `h`: specific enthalpy of air \[J/kg\]
    #[pyo3(name = "get_te_from_h")]
    #[staticmethod]
    pub fn get_te_from_h_py(h: f64) -> anyhow::Result<f64> {
        Ok(Self::get_te_from_h(h * uc::J_PER_KG)?.get::<si::degree_celsius>())
    }

)]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
pub struct Air {}
impl Init for Air {}
impl SerdeAPI for Air {}

impl Air {
    /// Returns density of air
    /// Source: <https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmet.html>  
    /// Note that if `None` is passed for either argument, function evaluation should be faster
    ///
    /// # Equations used
    /// - T = 15.04 - 0.00649 * h  
    /// - p = 101.29 * ((T + 273.1) / 288.08) ^ 5.256  
    ///
    /// # Arguments  
    /// * `te_air` - ambient temperature of air, defaults to 22 C
    /// * `h` - elevation above sea level, defaults to 180 m
    pub fn get_density(
        te_air: Option<si::ThermodynamicTemperature>,
        h: Option<si::Length>,
    ) -> si::MassDensity {
        let std_pressure_at_elev = |h: si::Length| -> si::Pressure {
            let std_temp_at_elev = (15.04 - 0.00649 * h.get::<si::meter>() + 273.15) * uc::KELVIN;
            (101.29e3 * uc::PASCAL)
                * ((std_temp_at_elev / (288.08 * uc::KELVIN))
                    .get::<si::ratio>()
                    .powf(5.256))
        };
        match (h, te_air) {
            (None, None) => *STD_DENSITY_AIR,
            (None, Some(te_air)) => *STD_PRESSURE_AIR / *R_AIR / te_air,
            (Some(h_val), None) => std_pressure_at_elev(h_val) / *R_AIR / *TE_STD_AIR,
            (Some(h_val), Some(te_air)) => std_pressure_at_elev(h_val) / *R_AIR / te_air,
        }
    }

    /// Returns thermal conductivity of air
    /// # Arguments
    /// - `te_air`: temperature of air
    pub fn get_therm_cond(
        te_air: si::ThermodynamicTemperature,
    ) -> anyhow::Result<si::ThermalConductivity> {
        Ok(
            THERMAL_CONDUCTIVITY_INTERP.interpolate(&[te_air.get::<si::kelvin>()])?
                * uc::WATT_PER_METER_KELVIN,
        )
    }

    /// Returns constant pressure specific heat of air
    /// # Arguments
    /// - `te_air`: temperature of air
    pub fn get_specific_heat_cp(
        te_air: si::ThermodynamicTemperature,
    ) -> anyhow::Result<si::SpecificHeatCapacity> {
        Ok(C_P_INTERP.interpolate(&[te_air.get::<si::kelvin>()])? * uc::J_PER_KG_K)
    }

    /// Returns specific enthalpy of air  
    /// # Arguments  
    /// - `te_air`: temperature of air
    pub fn get_specific_enthalpy(
        te_air: si::ThermodynamicTemperature,
    ) -> anyhow::Result<si::SpecificEnergy> {
        Ok(ENTHALPY_INTERP.interpolate(&[te_air.get::<si::kelvin>()])? * uc::J_PER_KG)
    }

    /// Returns thermal Prandtl number of air
    /// # Arguments
    /// - `te_air`: temperature of air     
    pub fn get_pr(te_air: si::ThermodynamicTemperature) -> anyhow::Result<si::Ratio> {
        Ok(PRANDTL_INTERP.interpolate(&[te_air.get::<si::kelvin>()])? * uc::R)
    }

    /// Returns dynamic viscosity \[Pa*s\] of air
    /// # Arguments
    /// te_air: temperature of air
    pub fn get_dyn_visc(
        te_air: si::ThermodynamicTemperature,
    ) -> anyhow::Result<si::DynamicViscosity> {
        Ok(DYN_VISC_INTERP.interpolate(&[te_air.get::<si::kelvin>()])? * uc::PASCAL_SECOND)
    }

    /// Returns temperature of air
    /// # Arguments
    /// `h`: specific enthalpy of air \[J/kg\]
    pub fn get_te_from_h(h: si::SpecificEnergy) -> anyhow::Result<si::ThermodynamicTemperature> {
        Ok(TEMP_FROM_ENTHALPY.interpolate(&[h.get::<si::joule_per_kilogram>()])? * uc::KELVIN)
    }
}

use air_static_props::*;

/// Fluid properties for calculations.  
///
/// Values obtained via (in Python):
/// ```python
/// >>> from CoolProp.CoolProp import PropsSI
/// >>> import numpy as np
/// >>> import pandas as pd
/// >>> T_degC = np.logspace(1, np.log10(5e3 + 70), 25) - 70
/// >>> T = T_degC + 273.15
/// >>> prop_dict = {
/// >>>     'T [°C]': T_degC,
/// >>>     'h [J/kg]': [0] * len(T),
/// >>>     'k [W/(m*K)]': [0] * len(T),
/// >>>     'rho [kg/m^3]': [0] * len(T),
/// >>>     'c_p [J/(kg*K)]': [0] * len(T),
/// >>>     'mu [Pa*s]': [0] * len(T),
/// >>> }
///
/// >>> for i, _ in enumerate(T_degC):
/// >>>     prop_dict['h [J/kg]'][i] = f"{PropsSI('H', 'P', 101325, 'T', T[i], 'Air'):.5g}" # specific enthalpy [J/(kg*K)]
/// >>>     prop_dict['k [W/(m*K)]'][i] = f"{PropsSI('L', 'P', 101325, 'T', T[i], 'Air'):.5g}" # thermal conductivity [W/(m*K)]
/// >>>     prop_dict['rho [kg/m^3]'][i] = f"{PropsSI('D', 'P', 101325, 'T', T[i], 'Air'):.5g}" # density [kg/m^3]
/// >>>     prop_dict['c_p [J/(kg*K)]'][i] = f"{PropsSI('C', 'P', 101325, 'T', T[i], 'Air'):.5g}" # density [kg/m^3]
/// >>>     prop_dict['mu [Pa*s]'][i] = f"{PropsSI('V', 'P', 101325, 'T', T[i], 'Air'):.5g}" # viscosity [Pa*s]
///
/// >>> prop_df = pd.DataFrame(data=prop_dict)
/// >>> pd.set_option('display.float_format', lambda x: '%.3g' % x)
/// >>> prop_df = prop_df.apply(np.float64)
/// ```
mod air_static_props {
    use super::*;
    lazy_static! {
        /// Array of temperatures at which properties are evaluated
        static ref TEMPERATURE_VALUES: Vec<si::ThermodynamicTemperature> = [-60.,
            -57.03690616,
            -53.1958198,
            -48.21658352,
            -41.7619528,
            -33.39475442,
            -22.54827664,
            -8.48788571,
            9.73873099,
            33.36606527,
            63.99440042,
            103.69819869,
            155.16660498,
            221.88558305,
            308.37402042,
            420.48979341,
            565.82652205,
            754.22788725,
            998.45434496,
            1315.04739396,
            1725.44993435,
            2257.45859876,
            2947.10642291,
            3841.10336915,
            5000.]
        .iter()
        .map(|x| *x * uc::KELVIN + *uc::CELSIUS_TO_KELVIN)
        .collect();
        pub static ref TEMP_FROM_ENTHALPY: Interp1D = Interp1D::new(
            ENTHALPY_VALUES.iter().map(|x| x.get::<si::joule_per_kilogram>()).collect::<Vec<f64>>(),
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
        /// Thermal conductivity values of air corresponding to temperature values
        static ref THERMAL_CONDUCTIVITY_VALUES: Vec<si::ThermalConductivity> = [
            0.019597,
            0.019841,
            0.020156,
            0.020561,
            0.021083,
            0.021753,
            0.022612,
            0.023708,
            0.025102,
            0.026867,
            0.02909,
            0.031875,
            0.035342,
            0.039633,
            0.044917,
            0.051398,
            0.059334,
            0.069059,
            0.081025,
            0.095855,
            0.11442,
            0.13797,
            0.16828,
            0.20795,
            0.26081]
        .iter()
        .map(|x| *x * uc::WATT_PER_METER_KELVIN)
        .collect();
        pub static ref THERMAL_CONDUCTIVITY_INTERP: Interp1D = Interp1D::new(
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            THERMAL_CONDUCTIVITY_VALUES.iter().map(|x| x.get::<si::watt_per_meter_degree_celsius>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
        /// Specific heat values of air corresponding to temperature values
        static ref C_P_VALUES: Vec<si::SpecificHeatCapacity> = [
            1006.2,
            1006.1,
            1006.,
            1005.9,
            1005.7,
            1005.6,
            1005.5,
            1005.6,
            1005.9,
            1006.6,
            1008.3,
            1011.6,
            1017.9,
            1028.9,
            1047.,
            1073.4,
            1107.6,
            1146.1,
            1184.5,
            1219.5,
            1250.1,
            1277.1,
            1301.7,
            1324.5,
            1347.]
        .iter()
        .map(|x| *x * uc::J_PER_KG_K)
        .collect();
        pub static ref C_P_INTERP: Interp1D = Interp1D::new(
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            C_P_VALUES.iter().map(|x| x.get::<si::joule_per_kilogram_kelvin>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
        static ref ENTHALPY_VALUES: Vec<si::SpecificEnergy> = [
            338940.,
            341930.,
            345790.,
            350800.,
            357290.,
            365710.,
            376610.,
            390750.,
            409080.,
            432860.,
            463710.,
            503800.,
            556020.,
            624280.,
            714030.,
            832880.,
            991400.,
            1203800.,
            1488700.,
            1869600.,
            2376700.,
            3049400.,
            3939100.,
            5113600.,
            6662000.]
        .iter()
        .map(|x| *x * uc::J_PER_KG)
        .collect();
        pub static ref ENTHALPY_INTERP: Interp1D = Interp1D::new(
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            ENTHALPY_VALUES.iter().map(|x| x.get::<si::joule_per_kilogram>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
        static ref DYN_VISCOSITY_VALUES: Vec<si::DynamicViscosity> = [
            1.4067e-05,
            1.4230e-05,
            1.4440e-05,
            1.4711e-05,
            1.5058e-05,
            1.5502e-05,
            1.6069e-05,
            1.6791e-05,
            1.7703e-05,
            1.8850e-05,
            2.0283e-05,
            2.2058e-05,
            2.4240e-05,
            2.6899e-05,
            3.0112e-05,
            3.3966e-05,
            3.8567e-05,
            4.4049e-05,
            5.0595e-05,
            5.8464e-05,
            6.8036e-05,
            7.9878e-05,
            9.4840e-05,
            1.1423e-04,
            1.4006e-04]
        .iter()
        .map(|x| *x * uc::PASCAL_SECOND)
        .collect();
        pub static ref DYN_VISC_INTERP: Interp1D = Interp1D::new(
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            DYN_VISCOSITY_VALUES.iter().map(|x| x.get::<si::pascal_second>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
        static ref PRANDTL_VALUES: Vec<si::Ratio> = DYN_VISCOSITY_VALUES
            .iter()
            .zip(C_P_VALUES.iter())
            .zip(THERMAL_CONDUCTIVITY_VALUES.iter())
            .map(|((mu, c_p), k)| -> si::Ratio {*mu * *c_p / *k})
            .collect::<Vec<si::Ratio>>();
        pub static ref PRANDTL_INTERP: Interp1D = Interp1D::new(
            TEMPERATURE_VALUES.iter().map(|x| x.get::<si::kelvin>()).collect::<Vec<f64>>(),
            PRANDTL_VALUES.iter().map(|x| x.get::<si::ratio>()).collect::<Vec<f64>>(),
            Strategy::Linear,
            Extrapolate::Error
        ).unwrap();
    }
}
