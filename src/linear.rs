use crate::conv::Convergence;
use ndarray::Array1;
pub mod naive_gauss;

#[derive(Debug, Clone, Default)]
pub struct LinearResult {
    pub x: Array1<f64>,
    pub conv: Convergence,
}
