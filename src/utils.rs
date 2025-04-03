use std::{error, fmt, result};

#[derive(Debug)]
pub struct KlError {
    msg: &'static str,
}
impl KlError {
    pub(crate) fn new(msg: &'static str) -> KlError {
        KlError { msg }
    }
    pub(crate) fn err<T>(msg: &'static str) -> KlResult<T> {
        Err(KlError::new(msg))
    }
    pub fn get_msg(&self) -> &'static str {
        self.msg
    }
}
impl fmt::Display for KlError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.msg)
    }
}
impl error::Error for KlError {}
pub type KlResult<T> = result::Result<T, KlError>;

/// A struct to represent a fraction.
/// The first vector represents the numerator,
/// and the second vector represents the denominator.
/// Terms are in ascending powers.
#[derive(Debug, PartialEq, Clone)]
pub struct Fraction(pub Vec<f64>, pub Vec<f64>);
impl fmt::Display for Fraction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}) / ({})", format_poly(&self.0), format_poly(&self.1))
    }
}
fn format_poly(p: &[f64]) -> String {
    let mut result = String::new();
    for (i, &coeff) in p.iter().enumerate() {
        if i > 0 {
            result.push_str(" + ");
        }
        result.push_str(&format_poly_term(coeff, i));
    }
    result
}
fn format_poly_term(coeff: f64, power: usize) -> String {
    let precision = 1;
    if power == 0 {
        format!("{:.*e}", precision, coeff)
    } else if power == 1 {
        format!("{:.*e}x", precision, coeff)
    } else {
        format!("{:.*e}x^{}", precision, coeff, power)
    }
}
/// Macro to create a fraction.
/// # Example
/// ```
/// use klsyn::{Fraction, fraction};
/// let f = fraction!([1, 2], [3, 4]);
/// assert_eq!(f, Fraction(vec![1.0, 2.0], vec![3.0, 4.0]));
/// ```
#[macro_export]
macro_rules! fraction {
    ($num:expr, $den:expr) => {
        Fraction(
            $num.iter().map(|&x| x as f64).collect(),
            $den.iter().map(|&x| x as f64).collect(),
        )
    };
}
