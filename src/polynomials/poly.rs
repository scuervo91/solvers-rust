use crate::SolverError;
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
}
