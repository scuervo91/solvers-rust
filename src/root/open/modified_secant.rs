use super::secant::secant;
use crate::RootResult;
use crate::SolverError;

pub fn modified_secant<F>(
    f: F,
    x0: f64,
    epsilon: f64,
    tol: f64,
    max_iter: usize,
) -> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
{
    let x1: f64 = if x0 == 0.0 {
        x0 + epsilon
    } else {
        x0 * (1.0 + epsilon)
    };

    return secant(f, x0, x1, tol, max_iter);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ok_modified_secant() {
        let fx = |x: f64| (-x).exp() - x;

        let result = modified_secant(fx, 0.0, 1e-4, 1e-6, 100);
        println!("Result: {:?}", result);
        assert!(result.is_ok());
    }
}
