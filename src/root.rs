use crate::conv::Convergence;

#[derive(Debug, Clone, Default)]
pub struct RootResult {
    pub root: f64,
    pub conv: Convergence,
    pub bracket: Option<(f64, f64)>,
}

pub mod bracket;
