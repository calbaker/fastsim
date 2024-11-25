use super::*;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, IsVariant)]
pub enum PowertrainType {
    ConventionalVehicle(Box<ConventionalVehicle>),
    HybridElectricVehicle(Box<HybridElectricVehicle>),
    BatteryElectricVehicle(Box<BatteryElectricVehicle>),
    // TODO: add PHEV here
}

impl SerdeAPI for PowertrainType {}
impl Init for PowertrainType {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::ConventionalVehicle(conv) => conv.init(),
            Self::HybridElectricVehicle(hev) => hev.init(),
            Self::BatteryElectricVehicle(bev) => bev.init(),
        }
    }
}

impl Powertrain for PowertrainType {
    fn set_curr_pwr_prop_out_max(
        &mut self,
        pwr_aux: si::Power,
        dt: si::Time,
        veh_state: &VehicleState,
    ) -> anyhow::Result<()> {
        match self {
            Self::ConventionalVehicle(v) => v.set_curr_pwr_prop_out_max(pwr_aux, dt, veh_state),
            Self::HybridElectricVehicle(v) => v.set_curr_pwr_prop_out_max(pwr_aux, dt, veh_state),
            Self::BatteryElectricVehicle(v) => v.set_curr_pwr_prop_out_max(pwr_aux, dt, veh_state),
        }
    }

    fn solve(
        &mut self,
        pwr_out_req: si::Power,
        veh_state: &VehicleState,
        enabled: bool,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        match self {
            Self::ConventionalVehicle(v) => v.solve(pwr_out_req, veh_state, enabled, dt),
            Self::HybridElectricVehicle(v) => v.solve(pwr_out_req, veh_state, enabled, dt),
            Self::BatteryElectricVehicle(v) => v.solve(pwr_out_req, veh_state, enabled, dt),
        }
    }

    fn solve_thermal(
        &mut self,
        te_amb: si::TemperatureInterval,
        heat_demand: si::Power,
        veh_speed: si::Velocity,
        dt: si::Time,
    ) -> anyhow::Result<()> {
        match self {
            Self::ConventionalVehicle(v) => v.solve_thermal(te_amb, heat_demand, veh_speed, dt),
            Self::HybridElectricVehicle(v) => v.solve_thermal(te_amb, heat_demand, veh_speed, dt),
            Self::BatteryElectricVehicle(v) => v.solve_thermal(te_amb, heat_demand, veh_speed, dt),
        }
    }

    fn get_curr_pwr_prop_out_max(&self) -> anyhow::Result<(si::Power, si::Power)> {
        match self {
            Self::ConventionalVehicle(v) => v.get_curr_pwr_prop_out_max(),
            Self::HybridElectricVehicle(v) => v.get_curr_pwr_prop_out_max(),
            Self::BatteryElectricVehicle(v) => v.get_curr_pwr_prop_out_max(),
        }
    }

    fn pwr_regen(&self) -> si::Power {
        match self {
            Self::ConventionalVehicle(v) => v.pwr_regen(),
            Self::HybridElectricVehicle(v) => v.pwr_regen(),
            Self::BatteryElectricVehicle(v) => v.pwr_regen(),
        }
    }
}

impl SaveInterval for PowertrainType {
    fn save_interval(&self) -> anyhow::Result<Option<usize>> {
        match self {
            PowertrainType::ConventionalVehicle(v) => v.save_interval(),
            PowertrainType::HybridElectricVehicle(v) => v.save_interval(),
            PowertrainType::BatteryElectricVehicle(v) => v.save_interval(),
        }
    }
    fn set_save_interval(&mut self, save_interval: Option<usize>) -> anyhow::Result<()> {
        match self {
            PowertrainType::ConventionalVehicle(v) => v.set_save_interval(save_interval),
            PowertrainType::HybridElectricVehicle(v) => v.set_save_interval(save_interval),
            PowertrainType::BatteryElectricVehicle(v) => v.set_save_interval(save_interval),
        }
    }
}

impl PowertrainType {
    pub fn conv_mut(&mut self) -> Option<&mut ConventionalVehicle> {
        match self {
            Self::ConventionalVehicle(conv) => Some(conv),
            _ => None,
        }
    }

    pub fn hev_mut(&mut self) -> Option<&mut HybridElectricVehicle> {
        match self {
            Self::HybridElectricVehicle(hev) => Some(hev),
            _ => None,
        }
    }

    // pub fn phev_mut(&mut self) -> Option<&mut> {
    //     self.pt_type.phev()
    // }

    pub fn bev_mut(&mut self) -> Option<&mut BatteryElectricVehicle> {
        match self {
            Self::BatteryElectricVehicle(bev) => Some(bev),
            _ => None,
        }
    }

    pub fn conv(&self) -> Option<&ConventionalVehicle> {
        match self {
            Self::ConventionalVehicle(conv) => Some(conv),
            _ => None,
        }
    }

    pub fn hev(&self) -> Option<&HybridElectricVehicle> {
        match self {
            Self::HybridElectricVehicle(hev) => Some(hev),
            _ => None,
        }
    }

    // pub fn phev(&self) -> Option<&> {
    //     self.pt_type.phev()
    // }

    pub fn bev(&self) -> Option<&BatteryElectricVehicle> {
        match self {
            Self::BatteryElectricVehicle(bev) => Some(bev),
            _ => None,
        }
    }

    pub fn fc(&self) -> Option<&FuelConverter> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => Some(&conv.fc),
            PowertrainType::HybridElectricVehicle(hev) => Some(&hev.fc),
            PowertrainType::BatteryElectricVehicle(_) => None,
        }
    }

    pub fn fc_mut(&mut self) -> Option<&mut FuelConverter> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => Some(&mut conv.fc),
            PowertrainType::HybridElectricVehicle(hev) => Some(&mut hev.fc),
            PowertrainType::BatteryElectricVehicle(_) => None,
        }
    }

    pub fn set_fc(&mut self, fc: FuelConverter) -> anyhow::Result<()> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => {
                conv.fc = fc;
                Ok(())
            }
            PowertrainType::HybridElectricVehicle(hev) => {
                hev.fc = fc;
                Ok(())
            }
            PowertrainType::BatteryElectricVehicle(_) => bail!("BEL has no FuelConverter."),
        }
    }

    pub fn fs(&self) -> Option<&FuelStorage> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => Some(&conv.fs),
            PowertrainType::HybridElectricVehicle(hev) => Some(&hev.fs),
            PowertrainType::BatteryElectricVehicle(_) => None,
        }
    }

    pub fn fs_mut(&mut self) -> Option<&mut FuelStorage> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => Some(&mut conv.fs),
            PowertrainType::HybridElectricVehicle(hev) => Some(&mut hev.fs),
            PowertrainType::BatteryElectricVehicle(_) => None,
        }
    }

    pub fn set_fs(&mut self, fs: FuelStorage) -> anyhow::Result<()> {
        match self {
            PowertrainType::ConventionalVehicle(conv) => {
                conv.fs = fs;
                Ok(())
            }
            PowertrainType::HybridElectricVehicle(hev) => {
                hev.fs = fs;
                Ok(())
            }
            PowertrainType::BatteryElectricVehicle(_) => bail!("BEL has no FuelConverter."),
        }
    }

    pub fn res(&self) -> Option<&ReversibleEnergyStorage> {
        match self {
            PowertrainType::ConventionalVehicle(_) => None,
            PowertrainType::HybridElectricVehicle(hev) => Some(&hev.res),
            PowertrainType::BatteryElectricVehicle(bev) => Some(&bev.res),
        }
    }

    pub fn res_mut(&mut self) -> Option<&mut ReversibleEnergyStorage> {
        match self {
            PowertrainType::ConventionalVehicle(_) => None,
            PowertrainType::HybridElectricVehicle(hev) => Some(&mut hev.res),
            PowertrainType::BatteryElectricVehicle(bev) => Some(&mut bev.res),
        }
    }

    pub fn set_res(&mut self, res: ReversibleEnergyStorage) -> anyhow::Result<()> {
        match self {
            PowertrainType::ConventionalVehicle(_) => {
                bail!("Conventional has no ReversibleEnergyStorage.")
            }
            PowertrainType::HybridElectricVehicle(veh) => {
                veh.res = res;
                Ok(())
            }
            PowertrainType::BatteryElectricVehicle(veh) => {
                veh.res = res;
                Ok(())
            }
        }
    }

    pub fn em(&self) -> Option<&ElectricMachine> {
        match self {
            PowertrainType::ConventionalVehicle(_conv) => None,
            PowertrainType::HybridElectricVehicle(hev) => Some(&hev.em),
            PowertrainType::BatteryElectricVehicle(bev) => Some(&bev.em),
        }
    }

    pub fn em_mut(&mut self) -> Option<&mut ElectricMachine> {
        match self {
            PowertrainType::ConventionalVehicle(_conv) => None,
            PowertrainType::HybridElectricVehicle(hev) => Some(&mut hev.em),
            PowertrainType::BatteryElectricVehicle(bev) => Some(&mut bev.em),
        }
    }

    pub fn set_em(&mut self, em: ElectricMachine) -> anyhow::Result<()> {
        match self {
            PowertrainType::ConventionalVehicle(_conv) => {
                Err(anyhow!("ConventionalVehicle has no `ElectricMachine`"))
            }
            PowertrainType::HybridElectricVehicle(hev) => {
                hev.em = em;
                Ok(())
            }
            PowertrainType::BatteryElectricVehicle(bev) => {
                bev.em = em;
                Ok(())
            }
        }
    }
}

impl SaveState for PowertrainType {
    fn save_state(&mut self) {
        match self {
            Self::ConventionalVehicle(conv) => conv.save_state(),
            Self::HybridElectricVehicle(hev) => hev.save_state(),
            Self::BatteryElectricVehicle(bev) => bev.save_state(),
        }
    }
}

impl Step for PowertrainType {
    fn step(&mut self) {
        match self {
            Self::ConventionalVehicle(conv) => conv.step(),
            Self::HybridElectricVehicle(hev) => hev.step(),
            Self::BatteryElectricVehicle(bev) => bev.step(),
        }
    }
}

#[allow(clippy::to_string_trait_impl)]
impl std::string::ToString for PowertrainType {
    fn to_string(&self) -> String {
        match self {
            PowertrainType::ConventionalVehicle(_) => String::from("Conv"),
            PowertrainType::HybridElectricVehicle(_) => String::from("HEV"),
            PowertrainType::BatteryElectricVehicle(_) => String::from("BEV"),
        }
    }
}
