struct Polynomial {
    coeffs: Vec<f64>
}

impl Polynomial {
    pub fn new(coeffs: Vec<f64>) -> Self {
        return Self { coeffs }
    }
}