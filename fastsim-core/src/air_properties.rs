use super::imports::*;
use super::*;

lazy_static! {
    pub static ref TE_AIR_DEFAULT: si::ThermodynamicTemperature = (22. + 273.15) * uc::KELVIN;
    pub static ref STD_PRESSURE_DEFAULT: si::Pressure = 99_346.3 * uc::PASCAL;
    pub static ref STD_DENSITY: si::MassDensity = 1.172 * uc::KGPM3;
    pub static ref GAS_CONSTANT: si::SpecificHeatCapacity = 287.0 * uc::J_PER_KG_K;
}

/// Returns density of air
/// Source: <https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmet.html>  
/// Note that if `None` is passed for either argument, function evaluation should be faster
///
/// # Equations used
/// T = 15.04 - .00649 * h  
/// p = 101.29 * [(T + 273.1)/288.08]^5.256  
///
/// # Arguments  
/// * `te_air` - ambient temperature of air, defaults to 22 C
/// * `h` - elevation above sea level, defaults to 180 m
pub fn get_density_air(
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
        (None, None) => *STD_DENSITY,
        (None, Some(te_air)) => *STD_PRESSURE_DEFAULT / *GAS_CONSTANT / te_air,
        (Some(h_val), None) => std_pressure_at_elev(h_val) / *GAS_CONSTANT / *TE_AIR_DEFAULT,
        (Some(h_val), Some(te_air)) => std_pressure_at_elev(h_val) / *GAS_CONSTANT / te_air,
    }
}

#[cfg(feature = "pyo3")]
#[pyfunction(name = "get_density_air")]
/// Returns density of air [kg/m^3]
/// Source: <https://www.grc.nasa.gov/WWW/K-12/rocket/atmosmet.html>  
///
/// # Equations used
/// T = 15.04 - .00649 * h  
/// p = 101.29 * [(T + 273.1)/288.08]^5.256  
///
/// # Arguments  
/// * `te_air_deg_c` - optional ambient temperature [Celsius] of air, defaults to 22 C
/// * `h_m` - optional elevation [m] above sea level, defaults to 180 m
pub fn get_density_air_py(te_air_deg_c: Option<f64>, h_m: Option<f64>) -> f64 {
    get_density_air(
        te_air_deg_c.map(|te_air_deg_c| (te_air_deg_c + 273.15) * uc::KELVIN),
        h_m.map(|h_m| h_m * uc::M),
    )
    .get::<si::kilogram_per_cubic_meter>()
}
