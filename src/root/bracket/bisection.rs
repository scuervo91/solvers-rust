pub fn bisection<F>(f: F,x0:f64, x1:f64) -> Result<f64, String>
where
    F: Fn(f64) -> f64
{
    let f_1 = f(x0)
    let f_2 = f(x1)
    let prod = f_1 * f_2

    if prod >= 0 {
        Err()
    } 

    
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisection() {
        let result = bisection(|x| x - 5.0);
        assert!(result.is_ok());
    }
}



