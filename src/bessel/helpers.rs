use crate::consts::{ROOT5_DBL_EPS, SQRT_DLB_EPS};
use crate::gamma::Gamma;
use crate::result::{SpecFunCode, SpecFunResult};

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
    (30375.0 * tpow[3] - 369603.0 * tpow[5] + 765765.0 * tpow[7] - 425425.0 * tpow[9]) / 414720.0
}

#[inline]
fn debye_u4(tpow: &[f64]) -> f64 {
    (4465125.0 * tpow[4] - 94121676.0 * tpow[6] + 349922430.0 * tpow[8] - 446185740.0 * tpow[10]
        + 185910725.0 * tpow[12])
        / 39813120.0
}

#[inline]
fn debye_u5(tpow: &[f64]) -> f64 {
    (1519035525.0 * tpow[5] - 49286948607.0 * tpow[7] + 284499769554.0 * tpow[9]
        - 614135872350.0 * tpow[11]
        + 566098157625.0 * tpow[13]
        - 188699385875.0 * tpow[15])
        / 6688604160.0
}

/// These are of use in calculating the oscillating Bessel functions.
// cos(y - pi/4 + eps)
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
//
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
        return result;
    } else {
        result.val = 0.0;
        result.err = 0.0;
        result
    }
}

// nu -> Inf; uniform in x > 0  [Abramowitz+Stegun, 9.7.8]
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
        let sum = 1.0 - debye_u1(&tpow) / nu + debye_u2(&tpow) / (nu * nu) -
            debye_u3(&tpow) / (nu * nu * nu) +
            debye_u4(&tpow) / (nu * nu * nu * nu) -
            debye_u5(&tpow) / (nu * nu * nu * nu * nu);
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

/* Evaluate J_mu(x),J_{mu+1}(x) and Y_mu(x),Y_{mu+1}(x)  for |mu| < 1/2
 */
pub(super) fn besseljyv_restricted(mu: f64, x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>,
                                                        SpecFunResult<f64>, SpecFunResult<f64>) {
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
    } else {
        if x < 2.0 {
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
            let (p, q) = bessel_jy_steed_cf2(mu, x);
            let jprime_j_ratio = mu / x - j_ratio.val;
            let gamma = (p - jprime_j_ratio) / q;
            jmu.val = j_sgn * (2.0 / (std::f64::consts::PI * x) / (q + gamma * (p -
                jprime_j_ratio)))
                .sqrt();
            jmu.err = 4.0 * f64::EPSILON * (jmu.val).abs();
            jmup1.val = j_ratio * jmu.val;
            jmup1.err = (j_ratio).abs() * jmu.err;
            ymu.val = gamma * jmu.val;
            ymu.err = (gamma).abs() * jmu.err;
            ymup1.val = ymu.val * (mu / x - p - q / gamma);
            ymup1.err = ymu.err * (mu / x - p - q / gamma).abs() + 4.0 * f64::EPSILON * ymup1.val.abs();
        } else {
            /* Use asymptotics for large argument. */
            jmu = besseljv_asympx_e(mu, x);
            jmup1 = besseljv_asympx_e(mu + 1.0, x);
            ymu = besselyv_asympx_e(mu, x);
            ymup1 = besselyv_asympx_e(mu + 1.0, x);
        }
    }

    (jmu, jmup1, ymu, ymup1)
}

pub(super) fn besselj_cf1(nu: f64, x: f64) -> (SpecFunResult<f64>, f64) {
    let recur_big = GSL_SQRT_DBL_MAX;
    let recur_small = GSL_SQRT_DBL_MIN;
    let maxiter = 10000;
    let mut n = 1;
    let mut anm2 = 1.0;
    let mut bnm2 = 0.0;
    let mut anm1 = 0.0;
    let mut bnm1 = 1.0;
    let mut a1 = x / (2.0 * (nu + 1.0));
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
            anm2 /= recur_big;
        } else if aan.abs() < recur_small || (bbn).abs() < recur_small {
            aan /= recur_small;
            bbn /= recur_small;
            anm1 /= recur_small;
            bnm1 /= recur_small;
            anm2 /= recur_small;
            bnm2 /= recur_small;
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

pub(super) fn besseljy_steed_cf2(nu: f64, x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let max_iter = 10000;
    let small = 1.0e-100;

    let mut p = SpecFunResult::<f64>::default();
    let mut q = SpecFunResult::<f64>::default();

    let mut i = 1;

    let mut x_inv = 1.0 / x;
    let mut a = 0.25 - nu * nu;
    p.val = -0.5 * x_inv;
    q.val = 1.0;
    let mut br = 2.0 * x;
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
        a += 2 * (i - 1) as f64;
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
