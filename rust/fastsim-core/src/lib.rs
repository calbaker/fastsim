// This needs to be a square logo to avoid stretching, and can have transparency
#![doc(html_logo_url = "https://www.nrel.gov/transportation/assets/images/icon-fastsim.jpg")]
//! Documentation for the Rust backend of the Future Automotive Systems Technology Simulator (FASTSim).

//! # Overview
//! FASTSim provides a simple way to compare powertrains and estimate the impact of technology
//! improvements on light-, medium-, and heavy-duty vehicle efficiency, performance, cost, and battery life.  
//! More information here: <https://www.nrel.gov/transportation/fastsim.html>

//! # Installation
//! Currently, the Rust backend is only available through a Python API.  
//! For installation instructions, see: <https://github.nrel.gov/MBAP/fastsim/blob/rust-port/fastsim/docs/README.md>

//! # Python Examples
//! ```python
//! import fastsim
//!
//! ## Load drive cycle by name
//! cyc_py = fastsim.cycle.Cycle.from_file("udds")
//! cyc_rust = cyc_py.to_rust()
//!
//! ## Load vehicle using database vehicle ID number
//! vnum = 1  
//! veh_py = fastsim.vehicle.Vehicle.from_vehdb(vnum)
//! veh_rust = veh_py.to_rust()
//!
//! ## Simulate
//! sd = fastsim.RustSimDrive(cyc_rust, veh_rust)
//! sd.sim_drive()
//! ```

extern crate ndarray;

#[macro_use]
pub mod macros;
extern crate proc_macros;
pub mod cycle;
pub mod params;
pub mod pyo3imports;
pub mod simdrive;
pub mod simdrive_impl;
pub mod utils;
pub mod vehicle;
