use crate::cheb::ChebSeries;
use crate::consts::{ROOT5_DBL_EPS, SQRT_DLB_EPS};
use crate::gamma::Gamma;
use crate::result::{SpecFunCode, SpecFunResult};
use lazy_static::lazy_static;
use std::f64::consts::PI;

lazy_static! {
    // nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu])
    static ref G1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            -1.145_164_083_662_683,
            0.006_360_853_113_470_843,
            0.001_862_451_930_072_068_4,
            0.000_152_833_085_873_453_5,
            0.000_017_017_464_011_802_038,
            -6.459_750_292_334_725e-7,
            -5.181_984_843_251_938e-8,
            4.518_909_289_485_818e-10,
            3.243_322_737_102_087e-11,
            6.830_943_402_494_752e-13,
            2.835_350_275_517_210_3e-14,
            -7.988_390_576_932_36e-16,
            -3.372_667_730_077_195e-17,
            -3.658_633_480_921_052e-20,
        ],
        a: -1.0,
        b: 1.0,
    };

    // nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu])
    static ref G2_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            1.882_645_524_949_671_9,
            -0.077_490_658_396_167_52,
            -0.018_256_714_847_324_93,
            0.000_633_803_020_907_489_6,
            0.000_076_229_054_350_872_9,
            -9.550_164_756_172_044e-7,
            -8.892_726_810_788_635e-8,
            -1.952_133_477_231_961_4e-9,
            -9.400_305_273_588_516e-11,
            4.687_513_384_953_239e-12,
            2.265_853_574_692_576e-13,
            -1.172_550_969_848_801_5e-15,
            -7.044_133_820_024_522e-17,
            -2.437_787_831_010_769_6e-18,
            -7.522_524_321_825_39e-20,
        ],
        a: -1.0,
        b: 1.0,
    };
}

#[inline]
fn debye_u1(tpow: &[f64]) -> f64 {
    (3.0 * tpow[1] - 5.0 * tpow[3]) / 24.0
}

#[inline]
fn debye_u2(tpow: &[f64]) -> f64 {
    (81.0 * tpow[2] - 462.0 * tpow[4] + 385.0 * tpow[6]) / 1152.0
}

#[inline]
fn debye_u3(tpow: &[f64]) -> f64 {
    (30375.0 * tpow[3] - 369_603.0 * tpow[5] + 765_765.0 * tpow[7] - 425_425.0 * tpow[9])
        / 414_720.0
}

#[inline]
fn debye_u4(tpow: &[f64]) -> f64 {
    (4_465_125.0 * tpow[4] - 94_121_676.0 * tpow[6] + 349_922_430.0 * tpow[8]
        - 446_185_740.0 * tpow[10]
        + 185_910_725.0 * tpow[12])
        / 39_813_120.0
}

#[inline]
fn debye_u5(tpow: &[f64]) -> f64 {
    (1_519_035_525.0 * tpow[5] - 49_286_948_607.0 * tpow[7] + 284_499_769_554.0 * tpow[9]
        - 614_135_872_350.0 * tpow[11]
        + 566_098_157_625.0 * tpow[13]
        - 188_699_385_875.0 * tpow[15])
        / 6_688_604_160.0
}

/// These are of use in calculating the oscillating Bessel functions.
// cos(y - pi/4 + eps)
#[allow(dead_code)]
pub(super) fn bessel_cos_pi4_e(y: f64, eps: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let sy = y.sin();
    let cy = y.cos();
    let s = sy + cy;
    let d = sy - cy;
    let abs_sum = cy.abs() + sy.abs();

    let (seps, ceps) = if eps.abs() < ROOT5_DBL_EPS {
        let e2 = eps * eps;
        (
            eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0)),
            1.0 - e2 / 2.0 * (1.0 - e2 / 12.0),
        )
    } else {
        (eps.sin(), eps.cos())
    };

    result.val = (ceps * s - seps * d) / std::f64::consts::SQRT_2;
    result.err =
        2.0 * f64::EPSILON * (ceps.abs() + seps.abs()) * abs_sum / std::f64::consts::SQRT_2;

    // Try to account for error in evaluation of sin(y), cos(y).
    // This is a little sticky because we don't really know
    // how the library routines are doing their argument reduction.
    // However, we will make a reasonable guess.
    // FIXME ?
    result.err *= if y > 1.0 / f64::EPSILON {
        0.5 * y
    } else if y > 1.0 / SQRT_DLB_EPS {
        256.0 * y * SQRT_DLB_EPS
    } else {
        1.0
    };

    result
}

/// These are of use in calculating the oscillating Bessel functions.
// sin(y - pi/4 + eps)
#[allow(dead_code)]
pub(super) fn bessel_sin_pi4_e(y: f64, eps: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let sy = y.sin();
    let cy = y.cos();
    let s = sy + cy;
    let d = sy - cy;
    let abs_sum = cy.abs() + sy.abs();

    let (seps, ceps) = if eps.abs() < ROOT5_DBL_EPS {
        let e2 = eps * eps;
        (
            eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0)),
            1.0 - e2 / 2.0 * (1.0 - e2 / 12.0),
        )
    } else {
        (eps.sin(), eps.cos())
    };

    result.val = (ceps * d + seps * s) / std::f64::consts::SQRT_2;
    result.err =
        2.0 * f64::EPSILON * (ceps.abs() + seps.abs()) * abs_sum / std::f64::consts::SQRT_2;

    // Try to account for error in evaluation of sin(y), cos(y).
    // This is a little sticky because we don't really know
    // how the library routines are doing their argument reduction.
    // However, we will make a reasonable guess.
    // FIXME ?
    result.err *= if y > 1.0 / f64::EPSILON {
        0.5 * y
    } else if y > 1.0 / SQRT_DLB_EPS {
        256.0 * y * SQRT_DLB_EPS
    } else {
        1.0
    };

    result
}

pub(super) fn bessel_ij_taylor_e(
    nu: f64,
    x: f64,
    sign: i32,
    kmax: i32,
    threshold: f64,
) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::default();

    if nu < 0.0 || x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.issue_warning(
            "bessel_ij_taylor_e",
            &[nu, x, sign as f64, kmax as f64, threshold],
        );
        result
    } else if x == 0.0 {
        if nu == 0.0 {
            result.val = 1.0;
        } else {
            result.val = 0.0;
        }
        result
    } else {
        let mut prefactor = SpecFunResult::default();
        let mut sum = SpecFunResult::<f64>::default();
        if nu == 0.0 {
            prefactor.val = 1.0;
            prefactor.err = 0.0;
            prefactor.code = SpecFunCode::Success;
        } else if nu < (i32::MAX - 1) as f64 {
            // Separate the integer part and use
            // y^nu/gamma(nu+1) = y^N / N! y^f/(N+1)_f to control the error.
            let n = (nu + 0.5).floor() as i32;
            let f = nu - n as f64;
            let poch_factor = (n as f64 + 1.0).poch_e(f);
            let tc_factor = crate::gamma::mono::taylorcoeff_e(n, 0.5 * x);
            let p = (0.5 * x).powf(f);
            prefactor.val = tc_factor.val * p / poch_factor.val;
            prefactor.err = tc_factor.err * p / poch_factor.val;
            prefactor.err += prefactor.val.abs() / poch_factor.val * poch_factor.err;
            prefactor.err += 2.0 * f64::EPSILON * prefactor.val.abs();
        } else {
            let lg = crate::gamma::mono::lngamma_e(nu + 1.0);
            let term1 = nu * (0.5 * x).ln();
            let term2 = lg.val;
            let ln_pre = term1 - term2;
            let ln_pre_err = f64::EPSILON * (term1.abs() + term2.abs()) + lg.err;
            prefactor = crate::exp::core::exp_err_e(ln_pre, ln_pre_err);
        }

        // Evaluate the sum
        {
            let y = sign as f64 * 0.25 * x * x;
            let mut sumk: f64 = 1.0;
            let mut term: f64 = 1.0;
            let mut k: usize = 0;
            for _k in 1..(kmax as usize + 1) {
                k = _k;
                term *= y / ((nu + k as f64) * k as f64);
                sumk += term;
                if (term / sumk).abs() < threshold {
                    break;
                }
            }

            sum.val = sumk;
            sum.err = threshold * sumk.abs();

            sum.code = if k >= kmax as usize {
                SpecFunCode::MaxIterErr
            } else {
                SpecFunCode::Success
            };
        }
        crate::elementary::multiply_err_e(prefactor.val, prefactor.err, sum.val, sum.err)
    }
}

/// Compute the asymptotic approximation to J_{nu}(x) using Hankel's
/// asymptotic expansion.
pub(super) fn besseljv_asympx_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let mu = 4.0 * nu * nu;
    let chi = x - (0.5 * nu + 0.25) * std::f64::consts::PI;
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut k: usize = 0;
    let mut t: f64 = 1.0;

    while k < 1000 {
        let kd = k as f64;

        let term = (mu - (2.0 * kd - 1.0) * (2.0 * kd - 1.0)) / (kd * (8.0 * x));

        t *= if k == 0 { 1.0 } else { -term };
        let convp: bool = t.abs() < f64::EPSILON * p.abs();
        p += t;

        k += 1;

        t *= term;
        let convq: bool = t.abs() < f64::EPSILON * q.abs();
        q += t;

        // To preserve the consistency of the series we need to exit when
        // p and q have the same number of terms.
        if convp && convq && kd > nu / 2.0 {
            break;
        }

        k += 1
    }

    {
        let pre = (2.0 / (x * std::f64::consts::PI)).sqrt();
        let c = chi.cos();
        let s = chi.sin();
        result.val = pre * (c * p - s * q);
        result.err =
            pre * f64::EPSILON * ((c * p).abs() + (s * q).abs() + t.abs()) * (1.0 + x.abs());
    }

    result
}

// x >> nu*nu+1
#[allow(dead_code)]
pub(super) fn besselyv_asympx_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let alpha = x;
    let beta = -0.5 * nu * std::f64::consts::PI;
    let ampl = 0f64; //bessel_asymp_Mv_e(nu, x);
    let theta = 0f64; //bessel_asymp_thetav_corr_e(nu, x);
    let sin_alpha = alpha.sin();
    let cos_alpha = alpha.cos();
    let sin_chi = (beta + theta).sin();
    let cos_chi = (beta + theta).cos();
    let sin_term = sin_alpha * cos_chi + sin_chi * cos_alpha;
    let sin_term_mag = (sin_alpha * cos_chi).abs() + (sin_chi * cos_alpha).abs();
    result.val = ampl * sin_term;
    result.err = ampl.abs() * f64::EPSILON * sin_term_mag;
    result.err += result.val.abs() * 2.0 * f64::EPSILON;

    if alpha.abs() > 1.0 / f64::EPSILON {
        result.err *= 0.5 * alpha.abs();
    } else if alpha.abs() > 1.0 / crate::consts::SQRT_DLB_EPS {
        result.err *= 256.0 * alpha.abs() * crate::consts::SQRT_DLB_EPS;
    }
    result
}

// x >> nu*nu+1
#[allow(dead_code)]
pub(super) fn besseliv_scaled_asympx_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let mu = 4.0 * nu * nu;
    let mum1 = mu - 1.0;
    let mum9 = mu - 9.0;
    let pre = 1.0 / (2.0 * std::f64::consts::PI * x).sqrt();
    let r = mu / x;
    result.val = pre * (1.0 - mum1 / (8.0 * x) + mum1 * mum9 / (128.0 * x * x));
    result.err = 2.0 * f64::EPSILON * result.val.abs() + pre * (0.1 * r * r * r).abs();
    result
}

#[allow(dead_code)]
pub(super) fn besselkv_scaled_asympx_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let mu = 4.0 * nu * nu;
    let mum1 = mu - 1.0;
    let mum9 = mu - 9.0;
    let pre = (std::f64::consts::PI / (2.0 * x)).sqrt();
    let r = nu / x;
    result.val = pre * (1.0 + mum1 / (8.0 * x) + mum1 * mum9 / (128.0 * x * x));
    result.err = 2.0 * f64::EPSILON * result.val.abs() + pre * (0.1 * r * r * r).abs();
    result
}

// nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.7]
//
// error:
//   The error has the form u_N(t)/nu^N  where  0 <= t <= 1.
//   It is not hard to show that |u_N(t)| is small for such t.
//   We have N=6 here, and |u_6(t)| < 0.025, so the error is clearly
//   bounded by 0.025/nu^6. This gives the asymptotic bound on nu
//   seen below as nu ~ 100. For general MACH_EPS it will be
//                     nu > 0.5 / MACH_EPS^(1/6)
//   When t is small, the bound is even better because |u_N(t)| vanishes
//   as t->0. In fact u_N(t) ~ C t^N as t->0, with C ~= 0.1.
//   We write
//                     err_N <= min(0.025, C(1/(1+(x/nu)^2))^3) / nu^6
//   therefore
//                     min(0.29/nu^2, 0.5/(nu^2+x^2)) < MACH_EPS^{1/3}
//   and this is the general form.
//
// empirical error analysis, assuming 14 digit requirement:
//   choose   x > 50.000 nu   ==>  nu >   3
//   choose   x > 10.000 nu   ==>  nu >  15
//   choose   x >  2.000 nu   ==>  nu >  50
//   choose   x >  1.000 nu   ==>  nu >  75
//   choose   x >  0.500 nu   ==>  nu >  80
//   choose   x >  0.100 nu   ==>  nu >  83
//
// This makes sense. For x << nu, the error will be of the form u_N(1)/nu^N,
// since the polynomial term will be evaluated near t=1, so the bound
// on nu will become constant for small x. Furthermore, increasing x with
// nu fixed will decrease the error.
#[allow(dead_code)]
pub(super) fn besseliv_scaled_asymp_unif_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let z = x / nu;
    let root_term = z.hypot(1.0);
    let pre = 1.0 / (2.0 * std::f64::consts::PI * nu * root_term).sqrt();
    let eta = root_term + (z / (1.0 + root_term)).ln();
    let ex_arg = if z < 1.0 / crate::consts::ROOT3_DBL_EPS {
        nu * (-z + eta)
    } else {
        -0.5 * nu / z * (1.0 - 1.0 / (12.0 * z * z))
    };
    let ex_result = crate::exp::core::exp_e(ex_arg);
    if ex_result.code == SpecFunCode::Success {
        let t = 1.0 / root_term;
        let sum;
        let mut tpow: [f64; 16] = [0f64; 16];
        tpow[0] = 1.0;
        for i in 1..16 {
            tpow[i] = t * tpow[i - 1];
        }
        sum = 1.0
            + debye_u1(&tpow) / nu
            + debye_u2(&tpow) / (nu * nu)
            + debye_u3(&tpow) / (nu * nu * nu)
            + debye_u4(&tpow) / (nu * nu * nu * nu)
            + debye_u5(&tpow) / (nu * nu * nu * nu * nu);
        result.val = pre * ex_result.val * sum;
        result.err = pre * ex_result.val / (nu * nu * nu * nu * nu * nu);
        result.err += pre * ex_result.err * sum.abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        result.val = 0.0;
        result.err = 0.0;
        result
    }
}

// nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
#[allow(dead_code)]
pub(super) fn besselkv_scaled_asymp_unif_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let z = x / nu;
    let root_term = z.hypot(1.0);
    let pre = (std::f64::consts::PI / (2.0 * nu * root_term)).sqrt();
    let eta = root_term + (z / (1.0 + root_term)).ln();
    let ex_arg = if z < 1.0 / crate::consts::ROOT3_DBL_EPS {
        nu * (z - eta)
    } else {
        0.5 * nu / z * (1.0 + 1.0 / (12.0 * z * z))
    };

    let ex_result = crate::exp::core::exp_e(ex_arg);
    if ex_result.code == SpecFunCode::Success {
        let t = 1.0 / root_term;
        let mut tpow = [0f64; 16];
        tpow[0] = 1.0;
        for i in 1..16 {
            tpow[i] = t * tpow[i - 1];
        }
        let sum = 1.0 - debye_u1(&tpow) / nu + debye_u2(&tpow) / (nu * nu)
            - debye_u3(&tpow) / (nu * nu * nu)
            + debye_u4(&tpow) / (nu * nu * nu * nu)
            - debye_u5(&tpow) / (nu * nu * nu * nu * nu);
        result.val = pre * ex_result.val * sum;
        result.err = pre * ex_result.err * sum.abs();
        result.err += pre * ex_result.val / (nu * nu * nu * nu * nu * nu);
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = 0.0;
        result.err = 0.0;
    }
    result
}

// Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
#[allow(dead_code)]
pub(super) fn besseljyv_restricted(
    mu: f64,
    x: f64,
) -> (
    SpecFunResult<f64>,
    SpecFunResult<f64>,
    SpecFunResult<f64>,
    SpecFunResult<f64>,
) {
    let mut jmu = SpecFunResult::<f64>::default();
    let mut jmup1 = SpecFunResult::<f64>::default();
    let mut ymu = SpecFunResult::<f64>::default();
    let mut ymup1 = SpecFunResult::<f64>::default();

    if x < 0.0 || mu.abs() > 0.5 {
        jmu.val = 0.0;
        jmu.err = 0.0;
        jmup1.val = 0.0;
        jmup1.err = 0.0;
        ymu.val = 0.0;
        ymu.err = 0.0;
        ymup1.val = 0.0;
        ymup1.err = 0.0;
    //GSL_ERROR("error", GSL_EDOM);
    } else if x == 0.0 {
        if mu == 0.0 {
            jmu.val = 1.0;
            jmu.err = 0.0;
        } else {
            jmu.val = 0.0;
            jmu.err = 0.0;
        }
        jmup1.val = 0.0;
        jmup1.err = 0.0;
        ymu.val = 0.0;
        ymu.err = 0.0;
        ymup1.val = 0.0;
        ymup1.err = 0.0;
    //GSL_ERROR("error", GSL_EDOM);
    } else if x < 2.0 {
        // Use Taylor series for J and the Temme series for Y.
        // The Taylor series for J requires nu > 0, so we shift
        // up one and use the recursion relation to get jmu, in
        // case mu < 0.
        jmup1 = bessel_ij_taylor_e(mu + 1.0, x, -1, 100, f64::EPSILON);
        let jmup2 = bessel_ij_taylor_e(mu + 2.0, x, -1, 100, f64::EPSILON);
        let c = 2.0 * (mu + 1.0) / x;
        jmu.val = c * jmup1.val - jmup2.val;
        jmu.err = c * jmup1.err + jmup2.err;
        jmu.err += 2.0 * f64::EPSILON * jmu.val.abs();
        let y_res = bessely_temme(mu, x);
        ymu = y_res.0;
        ymup1 = y_res.1;
    } else if x < 1000.0 {
        let (j_ratio, j_sgn) = besselj_cf1(mu, x);
        let (p, q) = besseljy_steed_cf2(mu, x);
        let jprime_j_ratio = mu / x - j_ratio.val;
        let gamma = (p.val - jprime_j_ratio) / q.val;
        jmu.val = j_sgn
            * (2.0 / (std::f64::consts::PI * x) / (q.val + gamma * (p.val - jprime_j_ratio)))
                .sqrt();
        jmu.err = 4.0 * f64::EPSILON * (jmu.val).abs();
        jmup1.val = j_ratio.val * jmu.val;
        jmup1.err = j_ratio.val.abs() * jmu.err;
        ymu.val = gamma * jmu.val;
        ymu.err = (gamma).abs() * jmu.err;
        ymup1.val = ymu.val * (mu / x - p.val - q.val / gamma);
        ymup1.err =
            ymu.err * (mu / x - p.val - q.val / gamma).abs() + 4.0 * f64::EPSILON * ymup1.val.abs();
    } else {
        /* Use asymptotics for large argument. */
        jmu = besseljv_asympx_e(mu, x);
        jmup1 = besseljv_asympx_e(mu + 1.0, x);
        ymu = besselyv_asympx_e(mu, x);
        ymup1 = besselyv_asympx_e(mu + 1.0, x);
    }

    (jmu, jmup1, ymu, ymup1)
}

#[allow(dead_code)]
pub(super) fn besselj_cf1(nu: f64, x: f64) -> (SpecFunResult<f64>, f64) {
    let recur_big = crate::consts::SQRT_DBL_MAX;
    let recur_small = crate::consts::SQRT_DBL_MIN;
    let maxiter = 10000;
    let mut n = 1;
    let mut anm2 = 1.0;
    let mut bnm2 = 0.0;
    let mut anm1 = 0.0;
    let mut bnm1 = 1.0;
    let a1 = x / (2.0 * (nu + 1.0));
    let mut aan = anm1 + a1 * anm2;
    let mut bbn = bnm1 + a1 * bnm2;
    let mut an;
    let mut ffn = aan / bbn;
    let mut dn = a1;
    let mut s = 1.0;

    while n < maxiter {
        n += 1;
        anm2 = anm1;
        bnm2 = bnm1;
        anm1 = aan;
        bnm1 = bbn;
        an = -x * x / (4.0 * (nu + n as f64 - 1.0) * (nu + n as f64));
        aan = anm1 + an * anm2;
        bbn = bnm1 + an * bnm2;

        if aan.abs() > recur_big || bbn.abs() > recur_big {
            aan /= recur_big;
            bbn /= recur_big;
            anm1 /= recur_big;
            bnm1 /= recur_big;
        } else if aan.abs() < recur_small || (bbn).abs() < recur_small {
            aan /= recur_small;
            bbn /= recur_small;
            anm1 /= recur_small;
            bnm1 /= recur_small;
        }

        let old_fn = ffn;
        ffn = aan / bbn;
        let del = old_fn / ffn;

        dn = 1.0 / (2.0 * (nu + n as f64) / x - dn);
        if dn < 0.0 {
            s = -s;
        }
        if (del - 1.0).abs() < 2.0 * f64::EPSILON {
            break;
        }
    }

    /* FIXME: we should return an error term here as well, because the
    error from this recurrence affects the overall error estimate. */

    let mut ratio = SpecFunResult::<f64>::default();
    ratio.val = ffn;
    let sgn = s;

    if n >= maxiter {
        ratio.code = SpecFunCode::MaxIterErr;
    }
    (ratio, sgn)
}

// Evaluate the continued fraction CF1 for I_{nu+1}/I_nu
// using Gautschi (Euler) equivalent series.
#[allow(dead_code)]
pub(super) fn besseli_cf1_ser(nu: f64, x: f64) -> SpecFunResult<f64> {
    let maxk = 20000;
    let mut tk = 1.0;
    let mut sum = 1.0;
    let mut rhok = 0.0;

    let mut result = SpecFunResult::<f64>::default();

    for k in 1..maxk {
        let ak = 0.25 * (x / (nu + k as f64)) * x / (nu + k as f64 + 1.0);
        rhok = -ak * (1.0 + rhok) / (1.0 + ak * (1.0 + rhok));
        tk *= rhok;
        sum += tk;
        if (tk / sum).abs() < f64::EPSILON {
            break;
        }
        if k == maxk - 1 {
            result.code = SpecFunCode::MaxIterErr;
        }
    }

    result.val = x / (2.0 * (nu + 1.0)) * sum;
    result
}

#[allow(dead_code)]
pub(super) fn besseljy_steed_cf2(nu: f64, x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let max_iter = 10000;
    let small = 1.0e-100;

    let mut p = SpecFunResult::<f64>::default();
    let mut q = SpecFunResult::<f64>::default();

    let mut i = 1;

    let x_inv = 1.0 / x;
    let mut a = 0.25 - nu * nu;
    p.val = -0.5 * x_inv;
    q.val = 1.0;
    let br = 2.0 * x;
    let mut bi = 2.0;
    let mut fact = a * x_inv / (p.val * p.val + q.val * q.val);
    let mut cr = br + q.val * fact;
    let mut ci = bi + p.val * fact;
    let mut den = br * br + bi * bi;
    let mut dr = br / den;
    let mut di = -bi / den;
    let mut dlr = cr * dr - ci * di;
    let mut dli = cr * di + ci * dr;
    let mut temp = p.val * dlr - q.val * dli;
    q.val = p.val * dli + q.val * dlr;
    p.val = temp;

    while i <= max_iter {
        i += 1;
        a += (2 * (i - 1)) as f64;
        bi += 2.0;
        dr = a * dr + br;
        di = a * di + bi;
        if dr.abs() + di.abs() < small {
            dr = small;
        }
        fact = a / (cr * cr + ci * ci);
        cr = br + cr * fact;
        ci = bi - ci * fact;
        if cr.abs() + ci.abs() < small {
            cr = small;
        }
        den = dr * dr + di * di;
        dr /= den;
        di /= -den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p.val * dlr - q.val * dli;
        q.val = p.val * dli + q.val * dlr;
        p.val = temp;
        if (dlr - 1.0).abs() + dli.abs() < f64::EPSILON {
            break;
        }
    }

    if i == max_iter {
        p.code = SpecFunCode::MaxIterErr;
        q.code = SpecFunCode::MaxIterErr;
    }
    (p, q)
}

/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 *
 * This is unstable for small x; x > 2 is a good cutoff.
 * Also requires |nu| < 1/2.
*/
#[allow(dead_code)]
pub(super) fn besselk_scaled_steed_temme_cf2(
    nu: f64,
    x: f64,
) -> (SpecFunResult<f64>, SpecFunResult<f64>, SpecFunResult<f64>) {
    let maxiter = 10000;

    let mut knu = SpecFunResult::<f64>::default();
    let mut knup1 = SpecFunResult::<f64>::default();
    let mut kpnu = SpecFunResult::<f64>::default();

    let mut bi = 2.0 * (1.0 + x);
    let mut di = 1.0 / bi;
    let mut delhi = di;
    let mut hi = di;

    let mut qi = 0.0;
    let mut qip1 = 1.0;

    let mut ai = -(0.25 - nu * nu);
    let a1 = ai;
    let mut ci = -ai;
    let mut qqi = -ai;

    let mut s = 1.0 + qqi * delhi;

    for i in 2..(maxiter + 1) {
        ai -= (2 * (i - 1)) as f64;
        ci = -ai * ci / i as f64;
        let tmp = (qi - bi * qip1) / ai;
        qi = qip1;
        qip1 = tmp;
        qqi += ci * qip1;
        bi += 2.0;
        di = 1.0 / (bi + ai * di);
        delhi *= bi * di - 1.0;
        hi += delhi;
        let dels = qqi * delhi;
        s += dels;
        if (dels / s).abs() < f64::EPSILON {
            break;
        }

        if i == maxiter - 1 {
            knu.code = SpecFunCode::MaxIterErr;
            knup1.code = SpecFunCode::MaxIterErr;
            kpnu.code = SpecFunCode::MaxIterErr;
        }
    }

    hi *= -a1;

    knu.val = (std::f64::consts::PI / (2.0 * x)).sqrt() / s;
    knup1.val = knu.val * (nu + x + 0.5 - hi) / x;
    kpnu.val = -knup1.val + nu / x * knu.val;

    (knu, knup1, kpnu)
}

pub(super) fn temme_gamma(nu: f64) -> (f64, f64, f64, f64) {
    let anu = nu.abs(); // functions are even
    let x: f64 = 4.0 * anu - 1.0;
    let rg1 = (*G1_CHEB).eval(x);
    let rg2 = (*G2_CHEB).eval(x);
    let g1 = rg1.val;
    let g2 = rg2.val;
    let g1mnu = 1.0 / (rg2.val + nu * rg1.val);
    let g1pnu = 1.0 / (rg2.val - nu * rg1.val);
    (g1pnu, g1mnu, g1, g2)
}

#[allow(dead_code)]
pub(super) fn bessely_temme(nu: f64, x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let mut ynu = SpecFunResult::<f64>::default();
    let mut ynup1 = SpecFunResult::<f64>::default();

    let max_iter = 15000;
    let half_x = 0.5 * x;
    let ln_half_x = half_x.ln();
    let half_x_nu = (nu * ln_half_x).exp();
    let pi_nu = std::f64::consts::PI * nu;
    let alpha = pi_nu / 2.0;
    let sigma = -nu * ln_half_x;
    let sinrat = if pi_nu.abs() < f64::EPSILON {
        1.0
    } else {
        pi_nu / pi_nu.sin()
    };
    let sinhrat = if sigma.abs() < f64::EPSILON {
        1.0
    } else {
        sigma.sinh() / sigma
    };
    let sinhalf = if alpha.abs() < f64::EPSILON {
        1.0
    } else {
        alpha.sin() / alpha
    };
    let sin_sqr = nu * PI * PI * 0.5 * sinhalf * sinhalf;

    let (g1pnu, g1mnu, g1, g2) = temme_gamma(nu);

    let mut fk = 2.0 / PI * sinrat * (sigma.cosh() * g1 - sinhrat * ln_half_x * g2);
    let mut pk = 1.0 / PI / half_x_nu * g1pnu;
    let mut qk = 1.0 / PI * half_x_nu * g1mnu;
    let mut hk;
    let mut ck = 1.0;

    let mut sum0 = fk + sin_sqr * qk;
    let mut sum1 = pk;

    let mut k = 0;

    while k < max_iter {
        k += 1;
        fk = (k as f64 * fk + pk + qk) / ((k * k) as f64 - nu * nu);
        ck *= -half_x * half_x / k as f64;
        pk /= k as f64 - nu;
        qk /= k as f64 + nu;
        let gk = fk + sin_sqr * qk;
        hk = -k as f64 * gk + pk;
        let del0 = ck * gk;
        let del1 = ck * hk;
        sum0 += del0;
        sum1 += del1;
        if del0.abs() < 0.5 * (1.0 + sum0.abs()) * f64::EPSILON {
            break;
        }
    }

    ynu.val = -sum0;
    ynu.err = (2.0 + 0.5 * k as f64) * f64::EPSILON * ynu.val.abs();
    ynup1.val = -sum1 * 2.0 / x;
    ynup1.err = (2.0 + 0.5 * k as f64) * f64::EPSILON * ynup1.val.abs();

    if k >= max_iter {
        ynu.code = SpecFunCode::MaxIterErr;
        ynup1.code = SpecFunCode::MaxIterErr;
    }

    (ynu, ynup1)
}

#[allow(dead_code)]
pub(super) fn besselk_scaled_temme(
    nu: f64,
    x: f64,
) -> (SpecFunResult<f64>, SpecFunResult<f64>, SpecFunResult<f64>) {
    let max_iter = 15000;

    let mut knu = SpecFunResult::<f64>::default();
    let mut knup1 = SpecFunResult::<f64>::default();
    let mut kpnu = SpecFunResult::<f64>::default();

    let half_x = 0.5 * x;
    let ln_half_x = half_x.ln();
    let half_x_nu = (nu * ln_half_x).exp();
    let pi_nu = PI * nu;
    let sigma = -nu * ln_half_x;
    let sinrat = if pi_nu.abs() < f64::EPSILON {
        1.0
    } else {
        pi_nu / pi_nu.sin()
    };
    let sinhrat = if sigma.abs() < f64::EPSILON {
        1.0
    } else {
        sigma.sinh() / sigma
    };
    let ex = x.exp();

    let mut k = 0;

    let (g1pnu, g1mnu, g1, g2) = temme_gamma(nu);

    let mut fk = sinrat * (sigma.cosh() * g1 - sinhrat * ln_half_x * g2);
    let mut pk = 0.5 / half_x_nu * g1pnu;
    let mut qk = 0.5 * half_x_nu * g1mnu;
    let mut hk = pk;
    let mut ck = 1.0;
    let mut sum0 = fk;
    let mut sum1 = hk;
    while k < max_iter {
        k += 1;
        fk = (k as f64 * fk + pk + qk) / ((k * k) as f64 - nu * nu);
        ck *= half_x * half_x / k as f64;
        pk /= k as f64 - nu;
        qk /= k as f64 + nu;
        hk = -k as f64 * fk + pk;
        let del0 = ck * fk;
        let del1 = ck * hk;
        sum0 += del0;
        sum1 += del1;
        if del0.abs() < 0.5 * sum0.abs() * f64::EPSILON {
            break;
        }
    }

    knu.val = sum0 * ex;
    knup1.val = sum1 * 2.0 / x * ex;
    kpnu.val = -knup1.val + nu / x * knu.val;

    if k == max_iter {
        knu.code = SpecFunCode::MaxIterErr;
        knup1.code = SpecFunCode::MaxIterErr;
        kpnu.code = SpecFunCode::MaxIterErr;
    }

    (knu, knup1, kpnu)
}

#[allow(dead_code)]
pub(super) fn besseljv_pos_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x == 0.0 {
        if nu == 0.0 {
            result.val = 1.0;
            result.err = 0.0;
        } else {
            result.val = 0.0;
            result.err = 0.0;
        }
    } else if x * x < 10.0 * (nu + 1.0) {
        return bessel_ij_taylor_e(nu, x, -1, 100, f64::EPSILON);
    } else if nu > 50.0 {
        return super::olver::besseljv_asymp_olver_e(nu, x);
    } else if x > 1000.0 {
        // We need this to avoid feeding large x to CF1; note that
        // due to the above check, we know that n <= 50. See similar
        // block in bessel_Jn.c.
        return besseljv_asympx_e(nu, x);
    } else {
        // -1/2 <= mu <= 1/2
        let nn = (nu + 0.5) as usize;
        let mu = nu - nn as f64;

        // Determine the J ratio at nu.
        let (jnup1_jnu, sgn_jnu) = besselj_cf1(nu, x);

        if x < 2.0 {
            // Determine Y_mu, Y_mup1 directly and recurse forward to nu.
            // Then use the CF1 information to solve for J_nu and J_nup1.
            let (ymu, ymup1) = bessely_temme(mu, x);

            let mut ynm1 = ymu.val;
            let mut yn = ymup1.val;
            let mut ynp1 = 0.0;

            for n in 1..nn {
                ynp1 = 2.0 * (mu + n as f64) / x * yn - ynm1;
                ynm1 = yn;
                yn = ynp1;
            }

            result.val = 2.0 / (std::f64::consts::PI * x) / (jnup1_jnu.val * yn - ynp1);
            result.err = f64::EPSILON * (nn as f64 + 2.0) * result.val.abs();
        } else {
            // Recurse backward from nu to mu, determining the J ratio
            // at mu. Use this together with a Steed method CF2 to
            // determine the actual J_mu, and thus obtain the normalization.

            let (p, q) = besseljy_steed_cf2(mu, x);

            let mut jnp1 = sgn_jnu * crate::consts::SQRT_DBL_MIN * jnup1_jnu.val;
            let mut jn = sgn_jnu * crate::consts::SQRT_DBL_MIN;

            for n in (1..(nn + 1)).rev() {
                let jnm1 = 2.0 * (mu + n as f64) / x * jn - jnp1;
                jnp1 = jn;
                jn = jnm1;
            }
            let jmup1_jmu = jnp1 / jn;
            let sgn_jmu = jn.signum();
            let jmuprime_jmu = mu / x - jmup1_jmu;

            let gamma = (p.val - jmuprime_jmu) / q.val;
            let jmu = sgn_jmu
                * (2.0 / (std::f64::consts::PI * x) / (q.val + gamma * (p.val - jmuprime_jmu)))
                    .sqrt();

            result.val = jmu * (sgn_jnu * crate::consts::SQRT_DBL_MIN) / jn;
            result.err = 2.0 * f64::EPSILON * (nn as f64 + 2.0) * result.val.abs();
        }
    }
    result
}

#[allow(dead_code)]
pub(super) fn besselyv_pos_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if nu > 50.0 {
        return super::olver::besselyv_asymp_olver_e(nu, x);
    } else {
        /* -1/2 <= mu <= 1/2 */
        let nn = (nu + 0.5) as usize;
        let mu = nu - nn as f64;

        let (ymu, ymup1) = if x < 2.0 {
            // Determine Ymu, Ymup1 directly. This is really
            // an optimization since this case could as well
            // be handled by a call to gsl_sf_bessel_JY_mu_restricted(),
            // as below.
            bessely_temme(mu, x)
        } else {
            // Determine Ymu, Ymup1 and Jmu, Jmup1.
            // &J_mu, &J_mup1, &Y_mu, &Y_mup1
            let res = besseljyv_restricted(mu, x);
            (res.2, res.3)
        };

        /* Forward recursion to get Ynu, Ynup1.
         */
        let mut ynm1 = ymu.val;
        let mut yn = ymup1.val;
        for n in 1..(nn + 1) {
            let ynp1 = 2.0 * (mu + n as f64) / x * yn - ynm1;
            ynm1 = yn;
            yn = ynp1;
        }

        result.val = ynm1;
        result.err = (nn as f64 + 1.0)
            * ynm1.abs()
            * ((ymu.err / ymu.val).abs() + (ymup1.err / ymup1.val).abs());
        result.err += 2.0 * f64::EPSILON * ynm1.abs();
    }
    result
}
