use crate::Convergence;
use crate::RootResult;
use crate::SolverError;

pub fn bisection<F>(
    f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<RootResult, SolverError>
where
    F: Fn(f64) -> f64,
{
    if a >= b {
        return Err(SolverError::InvalidInput(
            "a must be less than b".to_string(),
        ));
    }

    let mut fa: f64 = f(a);
    let fb: f64 = f(b);

    if fa.signum() * fb.signum() >= 0.0 {
        return Err(SolverError::NoRootInBrackets);
    }

    for iter in 1..=max_iter {
        let c: f64 = a + (b - a) / 2.0;

        let fc: f64 = f(c);

        let approx_error: f64 = ((b - a) / (b + a)).abs();

        if fc == 0.0 || approx_error < tol {
            let conv = Convergence {
                iterations: iter,
                nfev: 2 + iter, // f(a), f(b), then one f(c) per iteration
                residual_norm: fc.abs(),
                approx_error: approx_error,
            };
            let result = RootResult {
                root: c,
                conv,
                bracket: Some((a, b)),
            };
            return Ok(result);
        }

        if fa.signum() == fc.signum() {
            a = c;
            fa = fc;
        } else {
            b = c;
        }
    }

    return Err(SolverError::NotConverged { max_iter });
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

        let a: f64 = 12.0;
        let b: f64 = 16.0;
        let tol: f64 = 1e-6;
        let max_iter: usize = 100;

        let func = |c: f64| newton_second_law(c, g, m, t, v);

        let result = bisection(func, a, b, tol, max_iter);
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
            Err(SolverError::NotConverged { max_iter: 1 })
        ));
    }
}
