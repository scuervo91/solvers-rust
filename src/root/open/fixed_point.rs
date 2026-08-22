use crate::Convergence;
use crate::RootResult;
use crate::SolverError;

pub fn fixed_point<F>(f: F, x0: f64, tol: f64, max_iter: usize) -> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
{
    let mut xr: f64 = x0;
    let mut approx_error: f64 = 0.0;
    let mut residual_norm: f64 = 0.0;
    for iter in 1..=max_iter {
        let x_old = xr;
        xr = f(x_old);
        residual_norm = (xr - x_old).abs();

        approx_error = if xr != 0.0 {
            ((xr - x_old) / xr).abs()
        } else {
            residual_norm
        };

        if approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: iter,
                residual_norm,
                approx_error: approx_error,
            };
            let result = RootResult {
                root: xr,
                conv,
                ..Default::default()
            };
            return Ok(result);
        }
    }
    return Err(SolverError::NotConverged {
        max_iter,
        result: RootResult {
            root: xr,
            conv: Convergence {
                iterations: max_iter,
                nfev: max_iter,
                residual_norm,
                approx_error,
            },
            ..Default::default()
        },
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ok_fixed_point() {
        let result = fixed_point(|x| -x.exp(), 0.0, 1e-6, 100);
        println!("Result: {:?}", result);
        assert!(result.is_ok());
    }
}
