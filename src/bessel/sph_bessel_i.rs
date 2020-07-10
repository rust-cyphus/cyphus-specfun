use crate::result::{SpecFunCode, SpecFunResult};

#[allow(dead_code)]
pub(crate) fn sph_besselil_cf1(l: usize, x: f64, threshold: f64) -> (f64, SpecFunCode) {
    let kmax = 2000;
    let mut tk = 1.0;
    let mut sum = 1.0;
    let mut rhok = 0.0;
    let mut code = SpecFunCode::Success;
    let mut k = 1;

    loop {
        let ak = (x / (2 * l + 1 + 2 * k) as f64) * (x / (2 * l + 3 + 2 * k) as f64);
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk *= rhok;
        sum += tk;
        if (tk / sum).abs() < threshold {
            break;
        } else if k >= kmax {
            code = SpecFunCode::MaxIterErr;
            break;
        } else {
            k += 1;
        }
    }

    (x / (2 * l + 3) as f64 * sum, code)
}

#[allow(dead_code)]
pub(crate) fn sph_besseli0_scaled_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        result.val = 1.0;
        result
    } else if ax < 0.2 {
        let eax: f64 = (-ax).exp();
        let y: f64 = ax * ax;
        let c1: f64 = 1.0 / 6.0;
        let c2: f64 = 1.0 / 120.0;
        let c3: f64 = 1.0 / 5040.0;
        let c4: f64 = 1.0 / 362880.0;
        let c5: f64 = 1.0 / 39916800.0;
        let sum = c5
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
        result.val = eax * sum;
        result.err = 2.0 * f64::EPSILON * result.val;
        result
    } else if ax < -0.5 * crate::consts::LN_DBL_EPS {
        result.val = (1.0 - (-2.0 * ax).exp()) / (2.0 * ax);
        result.err = 2.0 * f64::EPSILON * result.val;
        result
    } else {
        result.val = 1.0 / (2.0 * ax);
        result.err = 2.0 * f64::EPSILON * result.val;
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_besseli1_scaled_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        result
    } else if ax < 3.0 * f64::MIN_POSITIVE {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::OverflowErr;
        result
    } else if ax < 0.25 {
        let eax: f64 = (-ax).exp();
        let y: f64 = ax * ax;
        let c1: f64 = 1.0 / 10.0;
        let c2: f64 = 1.0 / 280.0;
        let c3: f64 = 1.0 / 15120.0;
        let c4: f64 = 1.0 / 1330560.0;
        let c5: f64 = 1.0 / 172972800.0;
        let sum = c5
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
        result.val = eax * x / 3.0 * sum;
        result.err = 2.0 * f64::EPSILON * result.val;
        result
    } else {
        let ex = (-2.0 * ax).exp();
        result.val = 0.5 * (ax * (1.0 + ex) - (1.0 - ex)) / (ax * ax);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        if x < 0.0 {
            result.val *= -1.0;
        }
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_besseli2_scaled_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        result
    } else if ax < 4.0 * crate::consts::SQRT_DBL_MIN {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::OverflowErr;
        result
    } else if ax < 0.25 {
        let y: f64 = ax * ax;
        let c1: f64 = 1.0 / 14.0;
        let c2: f64 = 1.0 / 504.0;
        let c3: f64 = 1.0 / 33264.0;
        let c4: f64 = 1.0 / 3459456.0;
        let c5: f64 = 1.0 / 518918400.0;
        let sum = c5
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);
        let pre = (-ax).exp() * x * x / 15.0;
        result.val = pre * sum;
        result.err = 2.0 * f64::EPSILON * result.val;
        result
    } else {
        let ex = (-2.0 * ax).exp();
        let x2 = x * x;
        result.val = 0.5 * ((3.0 + x2) * (1.0 - ex) - 3.0 * ax * (1.0 + ex)) / (ax * ax * ax);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_besselil_scaled_e(l: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let mut sgn = 1.0;
    let ax = x.abs();
    let mut xx = x;

    if x < 0.0 {
        sgn = if l % 2 == 1 { -1.0 } else { 1.0 };
        xx *= -1.0;
    }

    if x == 0.0 {
        result.val = if l == 0 { 1.0 } else { 0.0 };
        result.err = 0.0;
        result
    } else if l == 0 {
        result = sph_besseli0_scaled_e(x);
        result.val *= sgn;
        result
    } else if l == 1 {
        result = sph_besseli1_scaled_e(x);
        result.val *= sgn;
        result
    } else if l == 2 {
        result = sph_besseli2_scaled_e(x);
        result.val *= sgn;
        result
    } else if x * x < 10.0 * (l as f64 + 1.5) / std::f64::consts::E {
        let b = super::bessel_helpers::bessel_ij_taylor_e(0.5 + l as f64, x, 1, 50, f64::EPSILON);
        let pre = (-ax).exp() * ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val = sgn * pre * b.val;
        result.err = pre * b.err;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else if l < 150 {
        let i0_scaled = sph_besseli0_scaled_e(ax);
        let (rat, _) = sph_besselil_cf1(l, ax, f64::EPSILON);
        let mut ilp1 = rat * crate::consts::SQRT_DBL_MIN;
        let mut il = crate::consts::SQRT_DBL_MIN;
        let mut ilm1;
        for ell in (1..=l).rev() {
            ilm1 = ilp1 + (2 * ell + 1) as f64 / x * il;
            ilp1 = il;
            il = ilm1;
        }
        result.val = sgn * i0_scaled.val * (crate::consts::SQRT_DBL_MIN / il);
        result.err = i0_scaled.err * (crate::consts::SQRT_DBL_MIN / il);
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else if (0.29 / (l * l + 1) as f64).min(0.5 / ((l * l + 1) as f64 + x * x))
        < 0.5 * crate::consts::ROOT3_DBL_MIN
    {
        result = super::bessel_helpers::besseliv_scaled_asymp_unif_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val *= sgn * pre;
        result.err *= pre;
        result
    } else {
        // recurse down from safe values
        let rt_term = ((0.5 * std::f64::consts::PI) / x).sqrt();
        let lmax = 2 + (1.2 / crate::consts::ROOT6_DBL_EPS) as usize;
        let r_ilp1 = super::bessel_helpers::besseliv_scaled_asymp_unif_e(0.5 + l as f64, x);
        let r_il = super::bessel_helpers::besseliv_scaled_asymp_unif_e(0.5 + l as f64, x);
        let mut ilp1 = r_ilp1.val;
        let mut il = r_il.val;
        let mut ilm1 = 0.0;
        ilp1 *= rt_term;
        il *= rt_term;
        for ell in ((l + 1)..=lmax).rev() {
            ilm1 = ilp1 + (2 * ell + 1) as f64 / x * il;
            ilp1 = il;
            il = ilm1;
        }
        result.val = sgn * ilm1;
        result.err = result.val.abs()
            * (f64::EPSILON + (r_ilp1.err / r_ilp1.val).abs() + (r_il.err / r_il.val).abs());
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}
