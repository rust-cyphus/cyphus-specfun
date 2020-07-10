use crate::result::{SpecFunCode, SpecFunResult};

#[allow(dead_code)]
pub(crate) fn sph_besselyl_small_x(l: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let den = crate::pow::powi_e(x, (l + 1) as i32);
    let num_fact = crate::gamma::mono::doublefact_e(2 * l - 1);

    if num_fact.code != SpecFunCode::Success || den.val == 0.0 {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result
    } else {
        let lmax = 200;
        let t: f64 = -0.5 * x * x;
        let mut sum: f64 = 1.0;
        let mut t_coeff: f64 = 1.0;
        let mut t_power: f64 = 1.0;
        let mut delta: f64;

        for i in 1..=lmax {
            t_coeff /= (i * (2 * (i - l) - 1)) as f64;
            t_power *= t;
            delta = t_power * t_coeff;
            sum += delta;
            if (delta / sum).abs() < 0.5 * f64::EPSILON {
                break;
            }
        }
        result.val = -num_fact.val / den.val * sum;
        result.err = f64::EPSILON * result.val.abs();
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_bessely0_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if 1.0 / std::f64::MAX > 0.0 && x < 1.0 / std::f64::MAX {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result
    } else {
        let cos = crate::trig::cos_e(x);
        result.val = -cos.val / x;
        result.err = (cos.err / x).abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_bessely1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if x < 1.0 / crate::consts::SQRT_DBL_MAX {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result
    } else if x < 0.25 {
        let y: f64 = x * x;
        let c1: f64 = 1.0 / 2.0;
        let c2: f64 = -1.0 / 8.0;
        let c3: f64 = 1.0 / 144.0;
        let c4: f64 = -1.0 / 5760.0;
        let c5: f64 = 1.0 / 403200.0;
        let c6: f64 = -1.0 / 43545600.0;
        let sum = c6
            .mul_add(y, c5)
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);

        result.val = -sum / y;
        result.err = f64::EPSILON * result.val.abs();
        result
    } else {
        let cos = crate::trig::cos_e(x);
        let sin = crate::trig::sin_e(x);
        result.val = -(cos.val / x + sin.val) / x;
        result.err = ((cos.err / x).abs() + sin.err) / x.abs();
        result.err += f64::EPSILON * ((sin.val / x).abs() + (cos.val / (x * x)).abs());
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_bessely2_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if x < 1.0 / crate::consts::SQRT_DBL_MAX {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result
    } else if x < 0.5 {
        let y: f64 = x * x;
        let c1: f64 = 1.0 / 6.0;
        let c2: f64 = 1.0 / 24.0;
        let c3: f64 = -1.0 / 144.0;
        let c4: f64 = 1.0 / 3456.0;
        let c5: f64 = -1.0 / 172800.0;
        let c6: f64 = 1.0 / 14515200.0;
        let c7: f64 = -1.0 / 1828915200.0;
        let sum = c7
            .mul_add(y, c6)
            .mul_add(y, c5)
            .mul_add(y, c4)
            .mul_add(y, c3)
            .mul_add(y, c2)
            .mul_add(y, c1)
            .mul_add(y, 1.0);

        result.val = -3.0 / (x * x * x) * sum;
        result.err = f64::EPSILON * result.val.abs();
        result
    } else {
        let cos = crate::trig::cos_e(x);
        let sin = crate::trig::sin_e(x);
        let a = 3.0 / (x * x);
        result.val = (1.0 - a) / x * cos.val - a * sin.val;
        result.err = cos.err * ((1.0 - a) / x).abs() + sin.err * a.abs();
        result.err += f64::EPSILON * ((cos.val / x).abs() + ((sin.val) / (x * x)).abs());
        result
    }
}

#[allow(dead_code)]
pub(crate) fn sph_besselyl_e(l: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if l == 0 {
        sph_bessely0_e(x)
    } else if l == 1 {
        sph_bessely1_e(x)
    } else if l == 2 {
        sph_bessely2_e(x)
    } else if x < 3.0 {
        sph_besselyl_small_x(l, x)
    } else if crate::consts::ROOT3_DBL_EPS * x > (l * l + l + 1) as f64 {
        result = super::bessel_helpers::besselyv_asympx_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val *= pre;
        result.err *= pre;
        result
    } else if l > 40 {
        result = super::olver::besselyv_asymp_olver_e(0.5 + l as f64, x);
        let pre = ((0.5 * std::f64::consts::PI) / x).sqrt();
        result.val *= pre;
        result.err *= pre;
        result
    } else {
        // recurse upward
        let r_by = sph_bessely1_e(x);
        let r_bym = sph_bessely0_e(x);
        let mut bym = r_bym.val;
        let mut by = r_by.val;
        let mut byp;
        for j in 1..l {
            byp = (2 * j + 1) as f64 / x * by - bym;
            bym = by;
            by = byp;
        }
        result.val = by;
        result.err = result.val.abs()
            * (f64::EPSILON + (r_by.err / r_by.val).abs() + (r_bym.err / r_bym.val).abs());
        result
    }
}
