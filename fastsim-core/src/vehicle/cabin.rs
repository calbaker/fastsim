use super::*;

/// Options for handling cabin thermal model
#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq, IsVariant)]
pub enum CabinOption {
    /// Basic single thermal capacitance cabin thermal model, including HVAC
    /// system and controls
    SingleCapacitanceCabin(Box<SingleCapacitanceCabin>),
    /// Cabin with interior and shell capacitances
    SingleCapacitanceCabinWithShell,
    /// no cabin thermal model
    #[default]
    None,
}
impl Init for CabinOption {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::SingleCapacitanceCabin(scc) => scc.init()?,
            Self::SingleCapacitanceCabinWithShell => {}
            Self::None => {}
        }
        Ok(())
    }
}
impl SerdeAPI for CabinOption {}

#[fastsim_api]
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, HistoryMethods)]
#[non_exhaustive]
/// Basic single thermal capacitance cabin thermal model, including HVAC
/// system and controls
pub struct SingleCapacitanceCabin {
    /// cabin shell thermal resistance \[m **2 * K / W\]
    pub cab_r_to_amb: f64,
    /// parameter for heat transfer coeff \[W / (m ** 2 * K)\] from cabin to ambient during
    /// vehicle stop
    pub cab_htc_to_amb_stop: f64,

    /// cabin thermal capacitance
    pub heat_capacitance: si::SpecificEnergy,
    /// cabin length, modeled as a flat plate
    pub length: si::Length,
    /// cabin width, modeled as a flat plate
    pub width: si::Length,
}

impl SerdeAPI for SingleCapacitanceCabin {}
impl Init for SingleCapacitanceCabin {}
