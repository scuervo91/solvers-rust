use thiserror::Error;

#[derive(Error, Debug)]
pub enum SolverError {
    // -- Exiting Errors --
    #[error("f(a) and f(b) must have opposite signs")]
    NoRootInBrackets,

    #[error("Invalid Input: {0}")]
    InvalidInput(String),

    #[error("Failed to converge within {max_iter} iterations")]
    NotConverged { max_iter: usize }

}