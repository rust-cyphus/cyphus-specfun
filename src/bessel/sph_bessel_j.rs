use crate::result::{SpecFunCode, SpecFunResult};

/// Compute the spherical bessel function of the first kind of order 0 along
/// with an error estimate.
#[allow(dead_code)]
pub(crate) fn sph_besselj0_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let ax = x.abs();

    if ax < 0.5 {
        let y = x * x;
        let c1: f64 = -1.0 / 6.0;
        let c2: f64 = 1.0 / 120.0;
        let c3: f64 = -1.0 / 5040.0;
        let c4: f64 = 1.0 / 362880.0;
        let c5: f64 = -1.0 / 39916800.0;
        let c6: f64 = 1.0 / 6227020800.0;
        result.val = c6
            .mul_add(y, c5)
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
    } else {
        result.val = x.sin() / x;
    }
    result.err = 2.0 * f64::EPSILON * result.val.abs();

    result
}

/// Compute the spherical bessel function of the first kind of order 1 along
/// with an error estimate.
#[allow(dead_code)]
pub(crate) fn sph_besselj1_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
    } else if ax < 3.1 * f64::MIN_POSITIVE {
        result.val = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else if ax < 0.25 {
        let y: f64 = x * x;
        let c1: f64 = -1.0 / 10.0;
        let c2: f64 = 1.0 / 280.0;
        let c3: f64 = -1.0 / 15120.0;
        let c4: f64 = 1.0 / 1330560.0;
        let c5: f64 = -1.0 / 172972800.0;
        let sum: f64 = c5
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
        result.val = x / 3.0 * sum;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        let cos = x.cos();
        let sin = x.sin();
        result.val = (sin / x - cos) / x;
        result.err = 2.0 * f64::EPSILON * ((sin / (x * x)).abs() + (cos / x).abs());
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the spherical bessel function of the first kind of order 2 along
/// with an error estimate.
#[allow(dead_code)]
pub(crate) fn sph_besselj2_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
    } else if ax < 4.0 * crate::consts::SQRT_DBL_MIN {
        result.val = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else if ax < 0.25 {
        let y: f64 = x * x;
        let c1: f64 = -1.0 / 14.0;
        let c2: f64 = 1.0 / 504.0;
        let c3: f64 = -1.0 / 33264.0;
        let c4: f64 = 1.0 / 3459456.0;
        let c5: f64 = -1.0 / 518918400.0;
        let c6: f64 = 1.0 / 105859353600.0;
        let c7: f64 = -1.0 / 28158588057600.0;
        let c8: f64 = 1.0 / 9461285587353600.0;
        let c9: f64 = -1.0 / 3916972233164390400.0;
        let sum: f64 = c9
            .mul_add(y, c8)
            .mul_add(y, c7)
            .mul_add(y, c6)
            .mul_add(y, c5)
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
        result.val = y / 15.0 * sum;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        let cos = x.cos();
        let sin = x.sin();
        let f = 3.0 / (x * x) - 1.0;
        result.val = (f * sin - 3.0 * cos / x) / x;
        result.err = 2.0 * f64::EPSILON * ((f * sin / x).abs() + 3.0 * (cos / (x * x)).abs());
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the spherical bessel function of the first kind of order `l` along
/// with an error estimate.
#[allow(dead_code)]
pub(crate) fn sph_besseljl_e(l: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x < 0.0 {
        result.val = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if l == 0 {
        sph_besselj0_e(x)
    } else if l == 1 {
        sph_besselj1_e(x)
    } else if l == 2 {
        sph_besselj2_e(x)
    } else if x * x < 10.0 * (0.5 + l as f64) / std::f64::consts::E {
        let b = super::bessel_helpers::bessel_ij_taylor_e(0.5 + l as f64, x, -1, 50, f64::EPSILON);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val = pre * b.val;
        result.err = pre * b.err;
        result.err += 2.0 * f64::EPSILON * result.val.abs();

        result
    } else if crate::consts::ROOT4_DBL_EPS * x > (l * l + l + 1) as f64 {
        let b = super::bessel_helpers::besseljv_asympx_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val = pre * b.val;
        result.err = 2.0 * f64::EPSILON * result.val + pre * b.err;
        result
    } else if l as f64 > 1.0 / crate::consts::ROOT6_DBL_EPS {
        let b = super::olver::besseljv_asymp_olver_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val = pre * b.val;
        result.err = 2.0 * f64::EPSILON * result.val + pre * b.err;
        result
    } else if x > 1000.0 && x > (l * l) as f64 {
        let b = super::bessel_helpers::besseljv_asympx_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val = pre * b.val;
        result.err = 2.0 * f64::EPSILON * result.val + pre * b.err;
        result
    } else {
        let (ratio, sgn) = super::bessel_helpers::besselj_cf1(0.5 + l as f64, x);
        let small = f64::MIN_POSITIVE / f64::EPSILON;
        let mut jlp1 = small * ratio.val;
        let mut jl = small;
        let mut jlm1;
        for ell in (1..(l + 1)).rev() {
            jlm1 = -jlp1 + (2 * ell + 1) as f64 / x * jl;
            jlp1 = jl;
            jl = jlm1;
        }

        if jl.abs() > jlp1.abs() {
            let j0 = sph_besselj0_e(x);
            let pre = small / jl;
            result.val = j0.val * pre;
            result.err = j0.err * pre.abs();
            result.err += 4.0 * f64::EPSILON * (0.5 * l as f64 + 1.0) * result.val.abs();
            result
        } else {
            let j1 = sph_besselj0_e(x);
            let pre = small / jlp1;
            result.val = j1.val * pre;
            result.err = j1.err * pre.abs();
            result.err += 4.0 * f64::EPSILON * (0.5 * l as f64 + 1.0) * result.val.abs();
            result
        }
    }
}
