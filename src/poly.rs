use crate::*;
use std::{cmp::max, vec};

/// Returns `true` if two polynomials are equal.
fn eq(p1: &[f32], p2: &[f32], eps: Option<f32>) -> bool {
    let eps = eps.unwrap_or_default().abs();

    let n1 = p1.len() - 1;
    let n2 = p2.len() - 1;
    let n = max(n1, n2);
    for i in 0..=n {
        let v1 = if i <= n1 { p1[i] } else { 0. };
        let v2 = if i <= n2 { p2[i] } else { 0. };
        if (v1 - v2).abs() > eps {
            return false;
        }
    }
    true
}

/// Adds two real polynomials.
fn add(p1: &[f32], p2: &[f32], eps: Option<f32>) -> KlResult<Vec<f32>> {
    let n1 = p1.len() - 1;
    let n2 = p2.len() - 1;
    let n3 = max(n1, n2);
    let mut result = vec![0.; n3 + 1];
    for i in 0..=n3 {
        let v1 = if i <= n1 { p1[i] } else { 0. };
        let v2 = if i <= n2 { p2[i] } else { 0. };
        result[i] = v1 + v2;
    }
    trim(&result, eps)
}

/// Multiplies two real polynomials.
fn multiply(p1: &[f32], p2: &[f32], eps: Option<f32>) -> KlResult<Vec<f32>> {
    if p1.is_empty() || p2.is_empty() {
        return KlError::err("zero length arrays");
    }
    if (p1.len() == 1 && p1[0] == 0.) || (p2.len() == 1 && p2[0] == 0.) {
        return Ok(vec![0.]);
    }
    let n1 = p1.len() - 1;
    let n2 = p2.len() - 1;
    let n3 = n1 + n2;
    let mut result = vec![0.; n3 + 1];
    for i in 0..=n3 {
        let mut t = 0.;
        let p1_start = if i > n2 { i - n2 } else { 0 };
        let p1_end = if i <= n1 { i } else { n1 };
        for j in p1_start..=p1_end {
            t += p1[j] * p2[i - j];
        }
        result[i] = t;
    }
    trim(&result, eps)
}

/// Divides two real polynomials.
/// Returns (quotient, remainder) = (p1 / p2, p1 % p2).
fn divide(p1: &[f32], p2: &[f32], eps: Option<f32>) -> KlResult<(Vec<f32>, Vec<f32>)> {
    if p1.is_empty() || p2.is_empty() {
        return KlError::err("zero length arrays");
    }
    let p1 = trim(p1, eps)?;
    let p2 = trim(p2, eps)?;
    if p2.len() == 1 {
        if p2[0] == 0. {
            return KlError::err("polynomial division by zero");
        }
        if p2[0] == 1. {
            return Ok((p1, vec![0.]));
        }
        let result = p1.iter().map(|&x| x / p2[0]).collect();
        return Ok((result, vec![0.]));
    }
    let n1 = p1.len() - 1;
    let n2 = p2.len() - 1;
    if n1 < n2 {
        return Ok((vec![0.], p1));
    }
    let mut a = p1;
    let lc2 = p2[n2]; // leading coefficient of p2
    for i in (0..=n1 - n2).rev() {
        let r = a[n2 + i] / lc2;
        a[n2 + i] = r;
        for j in 0..n2 {
            a[i + j] -= r * p2[j];
        }
    }
    let quotient = trim(&a[n2..], eps)?;
    let remainder = trim(&a[..n2], eps)?;
    Ok((quotient, remainder))
}

/// Returns the monic GCD (greatest common divisor) of two polynomials.
fn gcd(p1: &[f32], p2: &[f32], eps: Option<f32>) -> KlResult<Vec<f32>> {
    let mut r1 = trim(p1, eps)?;
    let mut r2 = trim(p2, eps)?;
    make_monic(&mut r1)?;
    make_monic(&mut r2)?;
    if r1.len() < r2.len() {
        std::mem::swap(&mut r1, &mut r2);
    }
    loop {
        if r2.len() < 2 {
            return Ok(vec![1.]); // GCD is 1
        }
        let mut remainder = divide(&r1, &r2, eps)?.1;
        if remainder.len() == 1 && remainder[0] == 0. {
            return Ok(r2);
        }
        make_monic(&mut remainder)?;
        r1 = r2;
        r2 = remainder;
    }
}

/// Trims top order zero coefficients.
fn trim(poly: &[f32], eps: Option<f32>) -> KlResult<Vec<f32>> {
    if poly.is_empty() {
        return KlError::err("zero length arrays");
    }
    let eps = eps.unwrap_or_default().abs();

    if poly.last().unwrap().abs() > eps {
        return Ok(poly.to_vec());
    }
    let len = poly.len();
    let mut new_len = len;
    while new_len > 0 && poly[new_len - 1].abs() <= eps {
        new_len -= 1;
    }
    if new_len == 0 {
        return Ok(vec![0.]);
    }
    Ok(poly[..new_len].to_vec())
}

/// Divides the coefficients by the leading coefficient.
fn make_monic(poly: &mut [f32]) -> KlResult<()> {
    let len = poly.len();
    if len == 0 {
        return KlError::err("zero length array");
    }
    let lc = poly[len - 1]; // leading coefficient
    if lc == 1. {
        return Ok(()); // already monic
    }
    if lc == 0. {
        return KlError::err("leading coefficient is zero"); // not trimmed
    }
    poly[len - 1] = 1.;
    poly.iter_mut().take(len - 1).for_each(|x| *x /= lc);
    Ok(())
}

/// Adds two fractions.
pub(crate) fn add_fractions(f1: &Fraction, f2: &Fraction, eps: Option<f32>) -> KlResult<Fraction> {
    let Fraction(n1, d1) = f1;
    let Fraction(n2, d2) = f2;

    if eq(d1, d2, eps) {
        return Ok(Fraction(add(n1, n2, eps)?, d1.to_vec()));
    }
    let g = gcd(d1, d2, eps)?;
    if g.len() == 1 && g[0] == 1. {
        let top = add(&multiply(n1, d2, eps)?, &multiply(n2, d1, eps)?, eps)?;
        let bottom = multiply(d1, d2, eps)?;
        return Ok(Fraction(top, bottom));
    }

    let (q1, _) = divide(d1, &g, eps)?;
    let (q2, _) = divide(d2, &g, eps)?;
    let top = add(&multiply(n1, &q2, eps)?, &multiply(n2, &q1, eps)?, eps)?;
    let bottom = multiply(d1, &q2, eps)?;
    Ok(Fraction(top, bottom))
}

/// Multiplies two fractions.
pub(crate) fn multiply_fractions(
    f1: &Fraction,
    f2: &Fraction,
    eps: Option<f32>,
) -> KlResult<Fraction> {
    let Fraction(n1, d1) = f1;
    let Fraction(n2, d2) = f2;
    let top = multiply(n1, n2, eps)?;
    let bottom = multiply(d1, d2, eps)?;
    Ok(Fraction(top, bottom))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fraction;

    #[test]
    fn test_eq() {
        assert!(eq(&[1., 2., 3.], &[1., 2., 3.], None));
        assert!(!eq(&[1., 2., 3.], &[1., 2., 4.], None));
        assert!(eq(&[1., 2., 3.], &[1., 2., 3., 0.], None));
        assert!(!eq(&[1., 2., 3.], &[1., 2., 3., 0.1], None));
    }

    #[test]
    fn test_add() {
        let result = add(&[1., 2., 3.], &[4., 5., 6.], None).unwrap();
        assert_eq!(result, vec![5., 7., 9.]);
    }

    #[test]
    fn test_multiply() {
        let result = multiply(&[1., 2.], &[3., 4.], None).unwrap();
        assert_eq!(result, vec![3., 10., 8.]);
    }

    #[test]
    fn test_divide() {
        let (quotient, remainder) = divide(&[2., 4., 6.], &[1., 2.], None).unwrap();
        assert_eq!(quotient, vec![0.5, 3.]);
        assert_eq!(remainder, vec![1.5]);
    }

    #[test]
    fn test_gcd() {
        let result = gcd(&[1., -3., 2.], &[1., -5.], None).unwrap();
        assert_eq!(result, vec![1.]);
        let result = gcd(&[1., -3., 2.], &[1., -3., 2.], None).unwrap();
        assert_eq!(result, vec![0.5, -1.5, 1.0]);
        let result = gcd(&[6., 5., -1.], &[6., 13., 8., 1.], None).unwrap();
        assert_eq!(result, vec![1., 1.]);
        let result = gcd(&[1.], &[0.], None);
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap().get_msg(),
            "leading coefficient is zero"
        );
    }

    #[test]
    fn test_trim() {
        let result = trim(&[1., 2., 3.], None).unwrap();
        assert_eq!(result, vec![1., 2., 3.]);
        let result = trim(&[1., 2., 0.], None).unwrap();
        assert_eq!(result, vec![1., 2.]);
        let result = trim(&[0., 0., 0.], None).unwrap();
        assert_eq!(result, vec![0.]);
        let result = trim(&[], None);
        assert!(result.is_err());
        assert_eq!(result.err().unwrap().get_msg(), "zero length arrays");
    }

    #[test]
    fn test_make_monic() {
        let mut poly = vec![6., 4., 2.];
        make_monic(&mut poly).unwrap();
        assert_eq!(poly, vec![3., 2., 1.]);
        let mut poly = vec![1., 0.];
        let result = make_monic(&mut poly);
        assert!(result.is_err());
        assert_eq!(
            result.err().unwrap().get_msg(),
            "leading coefficient is zero"
        );
    }

    #[test]
    fn test_fraction() {
        let f1 = fraction!([1, 2], [3, 4]);
        let f2 = fraction!([5., 6.], [7, 8]);
        let result = add_fractions(&f1, &f2, None).unwrap();
        assert_eq!(result, Fraction(vec![22., 60., 40.], vec![21., 52., 32.]));
        let result = multiply_fractions(&f1, &f2, None).unwrap();
        assert_eq!(result, Fraction(vec![5., 16., 12.], vec![21., 52., 32.]));
    }
}
