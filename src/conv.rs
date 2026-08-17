#[derive(Debug, Clone, Default)]
pub struct Convergence {
    pub iterations: usize,
    pub nfev: usize,
    pub residual_norm: f64,
    pub approx_error: f64,
}
