use crate::RootResult;
use crate::SolverError;
use crate::Convergence;


pub fn newton_raphson<F, DF>(
    f: F,
    df: DF,
    x0: f64,
    tol: f64,
    max_iter: usize
)-> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64
{
    let mut xr: f64 = x0;
    let mut approx_error: f64 = 0.0;
    let mut fx: f64 = 0.0;

    for iter in 1..=max_iter {
        let xr_old = xr;
        fx = f(xr_old);
        let dfx: f64 = df(xr_old);

        if dfx == 0.0 {
            return Err(SolverError::DivisionByZero);
        }

        xr = xr_old - fx / dfx;

        approx_error = if xr != 0.0 {
            ((xr - xr_old) / xr).abs()
        } else {
            (xr - xr_old).abs()
        };

        if approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: 2 * iter,
                residual_norm: fx.abs(),
                approx_error: approx_error
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
                nfev: 2 * max_iter,
                residual_norm: fx.abs(),
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
    fn test_ok_newton_raphson() {
        let fx = |x: f64| (-x).exp() - x;
        let dfx = |x: f64| -(-x).exp() - 1.0;

        let result = newton_raphson(fx, dfx, 0.0, 1e-6, 100);
        println!("Result: {:?}", result);
        assert!(result.is_ok());
    }
}



