use crate::Convergence;
use crate::RootResult;
use crate::SolverError;

pub fn bisection<F>(
    f: F,
    mut xl: f64,
    mut xu: f64,
    tol: f64,
    max_iter: usize,
) -> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
{
    if xl >= xu {
        return Err(SolverError::InvalidInput(
            "xl must be less than xu".to_string(),
        ));
    }

    let mut fl: f64 = f(xl);
    let fu: f64 = f(xu);

    if fl.signum() * fu.signum() >= 0.0 {
        return Err(SolverError::NoRootInBrackets);
    }

    let mut c: f64 = 0.0;
    let mut fc: f64 = 0.0;
    let mut approx_error: f64 = 0.0;

    for iter in 1..=max_iter {
        c = xl + (xu - xl) / 2.0;
        fc = f(c);
        approx_error = ((xu - xl) / (xu + xl)).abs();

        if fc == 0.0 || approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: 2 + iter, // f(xl), f(xu), then one f(c) per iteration
                residual_norm: fc.abs(),
                approx_error: approx_error,
            };
            let result = RootResult {
                root: c,
                conv,
                bracket: Some((xl, xu)),
            };
            return Ok(result);
        }

        if fl.signum() == fc.signum() {
            xl = c;
            fl = fc;
        } else {
            xu = c;
        }
    }

    let conv = Convergence {
        iterations: max_iter,
        nfev: 2 + max_iter,
        residual_norm: fc.abs(),
        approx_error,
    };
    let result = RootResult {
        root: c,
        conv,
        bracket: Some((xl, xu)),
    };
    return Err(SolverError::NotConverged { max_iter, result });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ok_bisection() {
        // Formula: f(c) = (g * m / c) * (1 - e^(-(c / m) * t)) - v
        // where:
        // g = gravitational constant
        // m = mass of the object
        // c = drag coefficient
        // t = time
        // v = velocity

        fn newton_second_law(c: f64, g: f64, m: f64, t: f64, v: f64) -> f64 {
            return (g * m / c) * (1.0 - (-(c / m) * t).exp()) - v;
        }
        let m: f64 = 68.1;
        let v: f64 = 40.0;
        let t: f64 = 10.0;
        let g: f64 = 9.81;

        let xl: f64 = 12.0;
        let xu: f64 = 16.0;
        let tol: f64 = 1e-6;
        let max_iter: usize = 100;

        let func = |c: f64| newton_second_law(c, g, m, t, v);

        let result = bisection(func, xl, xu, tol, max_iter);
        println!("Result: {:?}", result);
        assert!(result.is_ok());
    }

    #[test]
    fn test_bisection_a_greater_than_b() {
        let result = bisection(|x| x, 16.0, 12.0, 1e-6, 100);
        assert!(matches!(result, Err(SolverError::InvalidInput(_))));
    }

    #[test]
    fn test_bisection_no_root_in_brackets() {
        // f(x) = x^2 + 1 is positive on [0, 1]; no sign change
        let result = bisection(|x| x * x + 1.0, 0.0, 1.0, 1e-6, 100);
        assert!(matches!(result, Err(SolverError::NoRootInBrackets)));
    }

    #[test]
    fn test_bisection_not_converged() {
        // Valid bracket, but max_iter is too small to reach tol
        let result = bisection(|x| x - 1.0, 0.0, 3.0, 1e-12, 1);
        assert!(matches!(
            result,
            Err(SolverError::NotConverged { max_iter: 1, .. })
        ));
    }
}
