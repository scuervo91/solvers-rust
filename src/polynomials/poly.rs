use crate::RootResult;
use crate::SolverError;
use crate::conv::Convergence;
use std::cmp::max;
pub struct Polynomial {
    coeffs: Vec<f64>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<f64>) -> Self {
        return Self { coeffs };
    }

    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    pub fn order(&self) -> usize {
        return self.coeffs.len() - 1;
    }

    pub fn eval_at(&self, x: f64) -> f64 {
        return self
            .coeffs
            .iter()
            .fold(0.0, |acc, &coeff| (acc * x) + coeff);
    }

    pub fn eval(&self, x: &[f64]) -> Vec<f64> {
        return x.iter().map(|&x| self.eval_at(x)).collect();
    }

    pub fn eval_diff_at(&self, x: f64) -> (f64, f64) {
        return self.coeffs.iter().fold((0.0, 0.0), |(acc, dacc), &coeff| {
            ((acc * x) + coeff, (dacc * x) + acc)
        });
    }

    pub fn eval_diff(&self, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
        return x.iter().map(|&x| self.eval_diff_at(x)).unzip();
    }

    pub fn is_root(&self, x: f64, tol: Option<f64>) -> Result<bool, SolverError> {
        let tol = match tol {
            None => 1e-10,
            Some(t) if t.is_finite() && t > 0.0 => t,
            Some(_) => {
                return Err(SolverError::InvalidInput(
                    "Tolerance must be positive".to_string(),
                ));
            }
        };
        return Ok(self.eval_at(x).abs() <= tol);
    }

    pub fn deflate_monomial(&self, t: f64) -> Result<Polynomial, SolverError> {
        if self.coeffs.len() <= 1 {
            return Err(SolverError::InvalidInput(
                "Polynomial is too short to deflate".to_string(),
            ));
        }

        let n = self.coeffs.len();
        let is_root = self.is_root(t, None)?;

        if !is_root {
            return Err(SolverError::InvalidInput("t is not a root".to_string()));
        }

        // deflate by dividing by (x - t)

        let mut new_coeffs = Vec::with_capacity(self.coeffs.len().saturating_sub(1));

        let mut r = self.coeffs[0];

        for i in 1..n {
            let s = self.coeffs[i];
            let coeff_i = r;
            r = s + (r * t);
            new_coeffs.push(coeff_i);
        }

        return Ok(Polynomial::new(new_coeffs));
    }

    pub fn polsum(&self, other: &Polynomial) -> Result<Polynomial, SolverError> {
        let n = max(self.order(), other.order());

        let mut new_coeffs = Vec::with_capacity(n + 1);

        for i in 0..=n {
            let coeff_i = self.coeffs[i] + other.coeffs[i];
            new_coeffs.push(coeff_i);
        }

        return Ok(Polynomial::new(new_coeffs));
    }

    fn trim_leading_zeros(coeffs: &[f64]) -> Vec<f64> {
        match coeffs.iter().position(|&c| c != 0.0) {
            Some(i) => coeffs[i..].to_vec(),
            None => vec![0.0],
        }
    }

    pub fn poldiv(&self, other: &Polynomial) -> Result<(Polynomial, Polynomial), SolverError> {
        let dividend = Self::trim_leading_zeros(&self.coeffs);
        let divisor = Self::trim_leading_zeros(&other.coeffs);

        if divisor.len() == 1 && divisor[0] == 0.0 {
            return Err(SolverError::DivisionByZero);
        }

        let n = dividend.len() - 1;
        let m = divisor.len() - 1;

        if n < m {
            return Ok((Polynomial::new(vec![0.0]), Polynomial::new(dividend)));
        }

        let w = n - m;
        let mut r = dividend;
        let mut q = vec![0.0; w + 1];

        for i in 0..=w {
            q[i] = r[i] / divisor[0];
            for j in 1..=m {
                r[i + j] -= q[i] * divisor[j];
            }
        }

        let q = Self::trim_leading_zeros(&q);
        let rem = Self::trim_leading_zeros(&r[w + 1..]);
        return Ok((Polynomial::new(q), Polynomial::new(rem)));
    }

    // Root Finding Section
    pub fn muller_real_root(
        &self,
        mut x0: f64,
        mut x1: f64,
        mut x2: f64,
        tol: f64,
        max_iter: usize,
    ) -> Result<RootResult, SolverError> {
        let mut approx_error: f64 = 0.0;
        let mut x3: f64 = 0.0;

        for i in 1..=max_iter {
            let fx0 = self.eval_at(x0);
            let fx1 = self.eval_at(x1);
            let fx2 = self.eval_at(x2);

            let c = fx2;
            let h0 = x1 - x0;
            let h1 = x2 - x1;
            let delta0 = (fx1 - fx0) / h0;
            let delta1 = (fx2 - fx1) / h1;

            let a = (delta1 - delta0) / (h1 + h0);
            let b = (a * h1) + delta1;

            let disc = b * b - 4.0 * a * c;

            if disc < 0.0 {
                return Err(SolverError::InvalidInput(
                    "Discriminant is negative".to_string(),
                ));
            }

            x3 = x2 + ((-2.0 * c) / (b + b.signum() * disc.sqrt()));
            approx_error = ((x3 - x2) / x3).abs();

            if approx_error < tol {
                return Ok(RootResult {
                    root: x3,
                    conv: Convergence {
                        iterations: i,
                        nfev: i,
                        residual_norm: self.eval_at(x3).abs(),
                        approx_error: approx_error,
                    },
                    bracket: None,
                });
            }

            x0 = x1;
            x1 = x2;
            x2 = x3;
        }

        return Err(SolverError::NotConverged {
            max_iter: max_iter,
            result: RootResult {
                root: x3,
                conv: Convergence {
                    iterations: max_iter,
                    nfev: max_iter,
                    residual_norm: self.eval_at(x3).abs(),
                    approx_error: approx_error,
                },
                bracket: None,
            },
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // f5(x) = x^5 - 7x^4 - 3x^3 + 79x^2 - 46x - 120
    // coeffs are highest degree first (Horner)
    fn f5() -> Polynomial {
        return Polynomial::new(vec![1.0, -7.0, -3.0, 79.0, -46.0, -120.0]);
    }

    fn f2() -> Polynomial {
        Polynomial::new(vec![1.0, 2.0, -24.0])
    }

    fn f3() -> Polynomial {
        Polynomial::new(vec![1.0, 0.0, -13.0, -12.0])
    }

    #[test]
    fn test_eval_f5() {
        let p = f5();
        assert_eq!(p.order(), 5);
        assert_eq!(p.eval_at(0.0), -120.0);
        assert_eq!(p.eval_at(1.0), -96.0);
        assert_eq!(p.eval_at(2.0), 0.0); // root
    }

    #[test]
    fn test_eval() {
        let p = f5();
        let xs = [0.0, 1.0, 2.0];
        assert_eq!(p.eval(&xs), vec![-120.0, -96.0, 0.0]);
    }

    #[test]
    fn test_eval_diff_at_f5() {
        let p = f5();
        // f5'(x) = 5x^4 - 28x^3 - 9x^2 + 158x - 46
        assert_eq!(p.eval_diff_at(0.0), (-120.0, -46.0));
        assert_eq!(p.eval_diff_at(1.0), (-96.0, 80.0));
        assert_eq!(p.eval_diff_at(2.0), (0.0, 90.0));
    }

    #[test]
    fn test_eval_diff() {
        let p = f5();
        let xs = [0.0, 1.0, 2.0];
        let (ps, dps) = p.eval_diff(&xs);
        assert_eq!(ps, vec![-120.0, -96.0, 0.0]);
        assert_eq!(dps, vec![-46.0, 80.0, 90.0]);
    }

    #[test]
    fn test_deflate_monomial1() {
        let p = f2();

        let t: f64 = 4.0;

        let p_deflated = p.deflate_monomial(t).unwrap();
        assert_eq!(p_deflated.coeffs(), &[1.0, 6.0]);
        // (x - 4)(x + 6) = x^2 + 2x - 24
    }

    #[test]
    fn test_deflate_monomial2() {
        let p = f5();

        let t: f64 = -3.0;

        let p_deflated = p.deflate_monomial(t).unwrap();
        assert_eq!(p_deflated.coeffs(), &[1.0, -10.0, 27.0, -2.0, -40.0])
    }

    #[test]
    fn test_poldiv_exact() {
        let p = f2();
        let d = Polynomial::new(vec![1.0, -4.0]); // x - 4
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[1.0, 6.0]);
        assert_eq!(r.coeffs(), &[0.0]);
    }

    #[test]
    fn test_poldiv_remainder() {
        let p = f2();
        let d = Polynomial::new(vec![1.0, -1.0]); // x - 1
        let (q, r) = p.poldiv(&d).unwrap();
        // (x - 1)(x + 3) - 21 = x^2 + 2x - 24
        assert_eq!(q.coeffs(), &[1.0, 3.0]);
        assert_eq!(r.coeffs(), &[-21.0]);
    }

    #[test]
    fn test_poldiv_f5() {
        let p = f5();
        let d = Polynomial::new(vec![1.0, 3.0]); // x + 3
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[1.0, -10.0, 27.0, -2.0, -40.0]);
        assert_eq!(r.coeffs(), &[0.0]);
    }

    #[test]
    fn test_poldiv_lower_degree() {
        let p = f2();
        let d = f5();
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[0.0]);
        assert_eq!(r.coeffs(), p.coeffs());
    }

    #[test]
    fn test_poldiv_equal_degree() {
        let p = Polynomial::new(vec![2.0, 1.0]); // 2x + 1
        let d = Polynomial::new(vec![1.0, 1.0]); // x + 1
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[2.0]);
        assert_eq!(r.coeffs(), &[-1.0]);
    }

    #[test]
    fn test_poldiv_cubic() {
        let p = Polynomial::new(vec![1.0, -6.0, 11.0, -6.0]); // x^3 - 6x^2 + 11x - 6
        let d = Polynomial::new(vec![1.0, -1.0]); // x - 1
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[1.0, -5.0, 6.0]);
        assert_eq!(r.coeffs(), &[0.0]);
    }

    #[test]
    fn test_poldiv_remainder_trim() {
        // (x^3 - 1) / (x^2 + x + 1) = x - 1 remainder 0
        let p = Polynomial::new(vec![1.0, 0.0, 0.0, -1.0]);
        let d = Polynomial::new(vec![1.0, 1.0, 1.0]);
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[1.0, -1.0]);
        assert_eq!(r.coeffs(), &[0.0]);
        assert_eq!(r.order(), 0);
    }

    #[test]
    fn test_poldiv_leading_zeros() {
        let p = f2();
        let d = Polynomial::new(vec![0.0, 1.0, -4.0]); // 0x^2 + x - 4
        let (q, r) = p.poldiv(&d).unwrap();
        assert_eq!(q.coeffs(), &[1.0, 6.0]);
        assert_eq!(r.coeffs(), &[0.0]);
    }

    #[test]
    fn test_poldiv_zero_divisor() {
        let p = f2();
        let d = Polynomial::new(vec![0.0, 0.0]);
        let result = p.poldiv(&d);
        assert!(matches!(result, Err(SolverError::DivisionByZero)));
    }

    #[test]
    fn test_muller_real_root() {
        let p = f3();
        let result = p.muller_real_root(4.5, 5.5, 5.0, 1e-6, 100).unwrap();
        println!("result: {:?}", result);
        assert!((result.root - 4.0).abs() < 1e-6);
        assert!(result.conv.iterations < 100);
        assert!(result.conv.nfev < 100);
        assert!(result.conv.residual_norm < 1e-6);
        assert!(result.conv.approx_error < 1e-6);
        assert!(result.bracket.is_none());
    }
}
