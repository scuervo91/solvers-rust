pub struct Polynomial {
    coeffs: Vec<f64>,
}

impl Polynomial {
    pub fn new(coeffs: Vec<f64>) -> Self {
        return Self { coeffs };
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
}

#[cfg(test)]
mod tests {
    use super::*;

    // f5(x) = x^5 - 7x^4 - 3x^3 + 79x^2 - 46x - 120
    // coeffs are highest degree first (Horner)
    fn f5() -> Polynomial {
        Polynomial::new(vec![1.0, -7.0, -3.0, 79.0, -46.0, -120.0])
    }

    #[test]
    fn test_eval_f5() {
        let p = f5();
        assert_eq!(p.order(), 5);
        assert_eq!(p.eval_at(0.0), -120.0);
        assert_eq!(p.eval_at(1.0), -96.0);
        assert_eq!(p.eval_at(2.0), 0.0); // root
    }
}
