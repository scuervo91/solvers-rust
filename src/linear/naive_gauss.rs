use crate::SolverError;
use crate::conv::Convergence;
use crate::linear::LinearResult;
use ndarray::{Array1, Array2};

/// Solve `A x = b` by naive Gauss elimination (no pivoting).
///
/// 3×3 layout, textbook 1-based subscripts (`a_row,col`):
///
/// ```text
/// [ a11  a12  a13 ] [ x1 ]   [ b1 ]
/// [ a21  a22  a23 ] [ x2 ] = [ b2 ]
/// [ a31  a32  a33 ] [ x3 ]   [ b3 ]
/// ```
///
/// ndarray is 0-based: `a[[i, j]]` is `a_{i+1, j+1}`, `b[i]` is `b_{i+1}`.
///
/// ```text
/// a[[0, 0]] = a11    a[[0, 1]] = a12    a[[0, 2]] = a13    b[0] = b1
/// a[[1, 0]] = a21    a[[1, 1]] = a22    a[[1, 2]] = a23    b[1] = b2
/// a[[2, 0]] = a31    a[[2, 1]] = a32    a[[2, 2]] = a33    b[2] = b3
/// ```
pub fn naive_gauss(
    a_original: &Array2<f64>,
    b_original: &Array1<f64>,
) -> Result<LinearResult, SolverError> {
    let mut a = a_original.clone();
    let mut b = b_original.clone();

    // verify the shape of the matrix and vector are compatible
    if a.shape()[0] != a.shape()[1] {
        return Err(SolverError::InvalidShape(
            "Matrix must be square".to_string(),
        ));
    }
    if a.shape()[0] != b.len() {
        return Err(SolverError::InvalidShape(
            "Matrix and vector must have the same number of rows".to_string(),
        ));
    }

    let n: usize = a.shape()[0];

    let mut x = Array1::zeros(n);

    for k in 0..n {
        // forward elimination loop of the k-th column
        for i in k + 1..n {
            // compute factor for the i-th row
            let factor = a[[i, k]] / a[[k, k]];

            for j in k..n {
                a[[i, j]] -= factor * a[[k, j]]
            }
            b[i] -= factor * b[k]
        }
    }

    // backward substitution loop

    for k in (0..n).rev() {
        let mut sum_b = b[k];
        for i in k + 1..n {
            sum_b -= a[[k, i]] * x[i];
        }
        x[k] = sum_b / a[[k, k]];
    }

    let residual_norm = {
        let r = a_original.dot(&x) - b_original;
        r.iter().map(|v| v * v).sum::<f64>().sqrt()
    };

    return Ok(LinearResult {
        x,
        conv: Convergence {
            residual_norm,
            iterations: 1,
            nfev: 0,
            approx_error: 0.0,
        },
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    fn arr_1() -> (Array2<f64>, Array1<f64>) {
        //array 2d a
        let a = array![[3.0, -0.1, -0.2], [0.1, 7.0, -0.3], [0.3, -0.2, 10.0]];
        let b = array![7.85, -19.3, 71.4];
        (a, b)
    }

    #[test]
    fn test_ok_naive_gauss() {
        let (a, b) = arr_1();
        let result = naive_gauss(&a, &b).unwrap();
        // Chapra 9.5: x = [3, -2.5, 7]
        assert!((result.x[0] - 3.0).abs() < 1e-6);
        assert!((result.x[1] + 2.5).abs() < 1e-6);
        assert!((result.x[2] - 7.0).abs() < 1e-6);
        assert!(result.conv.residual_norm < 1e-10);
        assert_eq!(result.conv.iterations, 1);
        assert_eq!(result.conv.nfev, 0);
        assert_eq!(result.conv.approx_error, 0.0);
    }
}
