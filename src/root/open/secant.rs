use crate::Convergence;
use crate::RootResult;
use crate::SolverError;

pub fn secant<F>(
    f: F,
    x0: f64,
    x1: f64,
    tol: f64,
    max_iter: usize,
) -> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
{
    let mut xr_old: f64 = x0;
    let mut xr_new: f64 = x1;
    let mut xr: f64 = x1;

    let mut fx0: f64 = f(xr_old);
    let mut fx1: f64 = f(xr_new);
    let mut approx_error: f64 = 0.0;

    for iter in 1..=max_iter {
        if fx0 == fx1 {
            return Err(SolverError::DivisionByZero);
        }

        xr = xr_new - (fx1 * (xr_old - xr_new)) / (fx0 - fx1);
        let fx = f(xr);

        approx_error = if xr != 0.0 {
            ((xr - xr_new) / xr).abs()
        } else {
            (xr - xr_new).abs()
        };

        if approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: 2 + iter,
                residual_norm: fx.abs(),
                approx_error,
            };
            let result = RootResult {
                root: xr,
                conv,
                ..Default::default()
            };
            return Ok(result);
        }

        xr_old = xr_new;
        fx0 = fx1;
        xr_new = xr;
        fx1 = fx;
    }

    return Err(SolverError::NotConverged {
        max_iter,
        result: RootResult {
            root: xr,
            conv: Convergence {
                iterations: max_iter,
                nfev: 2 + max_iter,
                residual_norm: fx1.abs(),
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
    fn test_ok_secant() {
        let fx = |x: f64| (-x).exp() - x;

        let result = secant(fx, 0.0, 1.0, 1e-6, 100);
        println!("Result: {:?}", result);
        assert!(result.is_ok());
    }
}
