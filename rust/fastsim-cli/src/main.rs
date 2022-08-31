use clap::Parser;
use serde::{Deserialize, Serialize};

extern crate fastsim_core;
use fastsim_core::{cycle::RustCycle, simdrive::RustSimDrive, vehicle::RustVehicle};

/// Wrapper for fastsim.
/// After running `cargo build --release`, run with
/// ```bash
/// ./target/release/fastsim-cli --veh-file ~/Documents/GitHub/fastsim/fastsim/resources/vehdb/2012_Ford_Fusion.yaml --cyc-file ~/Documents/GitHub/fastsim/fastsim/resources/cycles/udds.csv
/// ```
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct FastSimApi {
    // `conflicts_with` tells rust one or the other, not both, can be provided
    /// Cycle as json string
    #[clap(long, value_parser, conflicts_with = "cyc_file")]
    cyc: Option<String>,
    #[clap(long, value_parser)]
    /// Path to cycle file (csv or yaml)
    cyc_file: Option<String>,
    #[clap(value_parser, long, conflicts_with = "veh_file")]
    /// Vehicle as json string
    veh: Option<String>,
    #[clap(long, value_parser)]
    /// Path to vehicle file (yaml)
    veh_file: Option<String>,
    #[clap(long, value_parser)]
    /// How to return results: `adopt_json`, `mpgge`, ... TBD
    res_fmt: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
struct AdoptResults {
    mpgge: f64,
    // add more results here
}

impl AdoptResults {
    pub fn to_json(&self) -> String {
        serde_json::to_string(&self).unwrap()
    }
}

pub fn main() {
    let fastsim_api = FastSimApi::parse();

    if let Some(_cyc_json_str) = fastsim_api.cyc {
        panic!("Need to implement: let cyc = RustCycle::from_json(cyc_json_str)");
    }

    let cyc = if let Some(cyc_file_path) = fastsim_api.cyc_file {
        RustCycle::from_file(&cyc_file_path)
    } else {
        RustCycle::from_file("../fastsim/resources/cycles/udds.csv")
    }
    .unwrap();

    let veh = RustVehicle::mock_vehicle();
    let mut sim_drive = RustSimDrive::new(cyc, veh);
    sim_drive.sim_drive(None, None).unwrap();

    let res_fmt = fastsim_api.res_fmt.unwrap_or_else(|| String::from("mpgge"));

    if res_fmt == "adopt_json" {
        let res = AdoptResults {
            mpgge: sim_drive.mpgge,
        };
        println!("{}", res.to_json());
    } else if res_fmt == "mpgge" {
        println!("{}", sim_drive.mpgge);
    } else {
        println!("Invalid option `{}` for `--res-fmt`", res_fmt);
    }
}

#[cfg(test)]
pub mod tests {
    use std::path::PathBuf;
    use std::process::Command;

    use assert_cmd::prelude::{CommandCargoExt, OutputAssertExt};
    use predicates::prelude::predicate;

    #[test]
    fn test_that_cli_app_produces_result() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("fastsim-cli")?;
        let mut cyc_file = project_root::get_project_root().unwrap();
        cyc_file.push(PathBuf::from("../fastsim/resources/cycles/udds.csv"));
        assert!(cyc_file.exists());
        let mut veh_file = project_root::get_project_root().unwrap();
        veh_file.push(PathBuf::from(
            "../fastsim/resources/vehdb/2012_Ford_Fusion.yaml",
        ));
        assert!(veh_file.exists());

        cmd.args([
            "--cyc-file",
            cyc_file.to_str().unwrap(),
            "--veh-file",
            veh_file.to_str().unwrap(),
        ]);
        cmd.assert()
            .success()
            .stdout(predicate::str::contains("32.4"));

        Ok(())
    }
}
