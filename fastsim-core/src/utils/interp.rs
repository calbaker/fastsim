use crate::imports::*;

/// Methods for proportionally scaling interpolator function data
pub trait InterpolatorMethods {
    fn set_min(&mut self, min: f64) -> anyhow::Result<()>;
    fn set_max(&mut self, max: f64) -> anyhow::Result<()>;
    fn set_range(&mut self, range: f64) -> anyhow::Result<()>;
}

impl InterpolatorMethods for Interpolator {
    #[allow(unused)]
    fn set_min(&mut self, min: f64) -> anyhow::Result<()> {
        let old_min = self.min()?;
        match self {
            Interpolator::Interp0D(value) => {
                *value = min;
                Ok(())
            }
            Interpolator::Interp1D(interp) => {
                todo!()
            }
            Interpolator::Interp2D(interp) => {
                todo!()
            }
            Interpolator::Interp3D(interp) => {
                todo!()
            }
            Interpolator::InterpND(interp) => {
                todo!()
            }
        }
    }

    fn set_max(&mut self, max: f64) -> anyhow::Result<()> {
        let old_max = self.max()?;
        match self {
            Interpolator::Interp0D(value) => {
                *value = max;
                Ok(())
            }
            Interpolator::Interp1D(interp) => {
                Ok(interp.set_f_x(interp.f_x().iter().map(|x| x * max / old_max).collect())?)
            }
            Interpolator::Interp2D(interp) => Ok(interp.set_f_xy(
                interp
                    .f_xy()
                    .iter()
                    .map(|v| v.iter().map(|x| x * max / old_max).collect())
                    .collect(),
            )?),
            Interpolator::Interp3D(interp) => Ok(interp.set_f_xyz(
                interp
                    .f_xyz()
                    .iter()
                    .map(|v0| {
                        v0.iter()
                            .map(|v1| v1.iter().map(|x| x * max / old_max).collect())
                            .collect()
                    })
                    .collect(),
            )?),
            Interpolator::InterpND(interp) => {
                Ok(interp.set_values(interp.values().map(|x| x * max / old_max))?)
            }
        }
    }

    fn set_range(&mut self, range: f64) -> anyhow::Result<()> {
        let old_max = self.max()?;
        let old_range = old_max - self.min()?;
        ensure!(old_range != 0., "Cannot modify range when min == max");
        match self {
            Interpolator::Interp0D(_value) => unreachable!("The above `ensure` should trigger"),
            Interpolator::Interp1D(interp) => Ok(interp.set_f_x(
                interp
                    .f_x()
                    .iter()
                    .map(|x| old_max + (x - old_max) * range / old_range)
                    .collect(),
            )?),
            Interpolator::Interp2D(interp) => Ok(interp.set_f_xy(
                interp
                    .f_xy()
                    .iter()
                    .map(|v| {
                        v.iter()
                            .map(|x| old_max + (x - old_max) * range / old_range)
                            .collect()
                    })
                    .collect(),
            )?),
            Interpolator::Interp3D(interp) => Ok(interp.set_f_xyz(
                interp
                    .f_xyz()
                    .iter()
                    .map(|v0| {
                        v0.iter()
                            .map(|v1| {
                                v1.iter()
                                    .map(|x| old_max + (x - old_max) * range / old_range)
                                    .collect()
                            })
                            .collect()
                    })
                    .collect(),
            )?),
            Interpolator::InterpND(interp) => Ok(interp.set_values(
                interp
                    .values()
                    .map(|x| old_max + (x - old_max) * range / old_range),
            )?),
        }
    }
}

impl Init for Interpolator {
    fn init(&mut self) -> anyhow::Result<()> {
        match self {
            Self::Interp0D(_) => {}
            Self::Interp1D(interp) => interp.validate()?,
            Self::Interp2D(interp) => interp.validate()?,
            Self::Interp3D(interp) => interp.validate()?,
            Self::InterpND(interp) => interp.validate()?,
        }
        Ok(())
    }
}
impl SerdeAPI for Interpolator {
    const RESOURCE_PREFIX: &'static str = "interpolators";
}
