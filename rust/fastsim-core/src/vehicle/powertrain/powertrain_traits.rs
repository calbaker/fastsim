//! Traits defining power flow interfaces for electric machines
use super::*;

// TODO: fix and uncomment or delete this
// pub trait ElectricMachine {
//     /// Sets current max power output given `pwr_in_max` from upstream component
//     fn set_cur_pwr_max_out(
//         &mut self,
//         pwr_in_max: si::Power,
//         pwr_aux: Option<si::Power>,
//     ) -> anyhow::Result<()>;
//     /// Sets current max power output rate given `pwr_rate_in_max` from upstream component
//     fn set_pwr_rate_out_max(&mut self, pwr_rate_in_max: si::PowerRate);
// }

pub trait Mass {
    /// Returns mass of Self, including contribution from any fields that implement `Mass`
    fn mass(&self) -> anyhow::Result<Option<si::Mass>>;

    /// Sets component mass to `mass`, or if `None` is provided for `mass`,
    /// sets mass based
    /// on other component parameters (e.g. power and power density)
    fn update_mass(&mut self, mass: Option<si::Mass>) -> anyhow::Result<()>;

    /// Checks if mass is consistent with other parameters
    fn check_mass_consistent(&self) -> anyhow::Result<()>;
}

/// Provides functions for solving powertrain
pub trait Powertrain {
    /// Returns maximum possible tractive power this component/system can produce, accounting for any aux power required.
    /// # Arguments
    /// - `pwr`
    /// - `dt`: time step size
    fn get_curr_pwr_out_max(
        &mut self,
        pwr_aux: si::Power,
        dt: si::Time,
    ) -> anyhow::Result<si::Power>;
    /// Solves for this powertrain system/component efficiency and sets/returns power output values.
    /// # Arguments
    /// - `pwr_out_req`: tractive power output required to achieve presribed speed  
    /// - `pwr_aux`: component-specific aux power demand (e.g. mechanical power if from engine/FC)  
    /// - `enabled`: whether component is actively running   
    /// - `dt`: time step size  
    fn solve(
        &mut self,
        pwr_out_req: si::Power,
        pwr_aux: si::Power,
        enabled: bool,
        dt: si::Time,
        assert_limits: bool,
    ) -> anyhow::Result<()>;
}
