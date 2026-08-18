use crate::RootResult;
use crate::SolverError;
use crate::Convergence;

pub fn fixed_point<F>(
    f: F,
    x0: f64,
    tol: f64,
    max_iter: usize,
) -> Result<RootResult, SolverError>
where 
    F: Fn(f64) -> f64,
{
    let mut xr: f64 = x0;
    let mut approx_error: f64 = 0.0;
    for iter in 1..=max_iter {
        let x_old = xr;
        xr = f(x_old);
        
        if xr != 0.0 {
            approx_error = ((xr - x_old)/xr).abs();
        }

        if approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: iter,
                residual_norm: xr.abs(),
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
    return Err(SolverError::NotConverged { max_iter, result: RootResult { root: xr, conv: Convergence { iterations: max_iter, nfev: max_iter, residual_norm: xr.abs(), approx_error: approx_error }, ..Default::default() } });
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