// Needed data from gamma mod
use crate::gamma::data::{BERN, GAMMA_XMAX};
use crate::result::{SpecFunCode, SpecFunResult};

// Need gamma functions
use crate::gamma::mono::gammainv_e;
use crate::gamma::mono::{lngamma_e, lngamma_sgn_e};
use crate::gamma::psi::digamma_e;

// Exponential functions
use crate::exp::core::{exp_err_e, expm1_e};
use crate::logarithm::ln_p1_e;

// Useful constants
use crate::consts::{LN_DBL_EPS, SQRT_3, SQRT_DBL_MIN};
use std::f64::consts::{E, LN_2, PI, SQRT_2};

// ((a)_x - 1)/x in the "small x" region where
// cancellation must be controlled.
//
// Based on SLATEC DPOCH1().
//
// When ABS(X) is so small that substantial cancellation will occur if
// the straightforward formula is used, we use an expansion due
// to Fields and discussed by Y. L. Luke, The Special Functions and Their
// Approximations, Vol. 1, Academic Press, 1969, page 34.
//
// The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
//        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
// In order to maintain significance in POCH1, we write for positive a
//        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
//                       = 1.0 + Q*EXPREL(Q) .
// Likewise the polynomial is written
//        POLY = 1.0 + X*POLY1(A,X) .
// Thus,
//        POCH1(A,X) = (POCH(A,X) - 1) / X
//                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
fn pochrel_smallx(a: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let sqrt_big = 1.0 / (2.0 * SQRT_2 * SQRT_3 * SQRT_DBL_MIN);
    let alneps = LN_DBL_EPS - LN_2;

    if x == 0.0 {
        digamma_e(a)
    } else {
        let bp = if a < -0.5 { 1.0 - a - x } else { a };
        let incr: i32 = if bp < 10.0 { (11.0 - bp) as i32 } else { 0 };
        let b = bp + incr as f64;

        let var = b + 0.5 * (x - 1.0);
        let alnvar = var.ln();
        let q = x * alnvar;

        let mut poly1 = 0.0;

        if var < sqrt_big {
            let nterms = (-0.5 * alneps / alnvar + 1.0) as usize;
            let var2 = 1.0 / var.powi(2);
            let rho = 0.5 * (x + 1.0);
            let mut term = var2;

            let mut gbern: [f64; 24] = [0.0; 24];
            gbern[1] = 1.0;
            gbern[2] = -rho / 12.0;
            poly1 = gbern[2] * term;

            if nterms > 20 {
                // nterms is too big, mayber precision is bad?
                result.code = SpecFunCode::SanityCheckErr;
                return result;
            }

            for k in 2..(nterms + 1) {
                let mut gbk = 0.0;
                for j in 1..(k + 1) {
                    gbk += BERN[k - j + 1] * gbern[j];
                }
                gbern[k + 1] = -rho * gbk / (k as f64);

                term *= ((2 * k) as f64 - 2.0 - x) * ((2 * k) as f64 - 1.0 - x) * var2;
                poly1 += gbern[k + 1] * term;
            }
        }

        let mut dexprl = expm1_e(q);
        if dexprl.code != SpecFunCode::Success {
            result.code = dexprl.code;
            return result;
        }
        dexprl.val /= q;
        poly1 *= x - 1.0;
        let mut poch1 = dexprl.val * (alnvar + q * poly1) + poly1;
        for i in (0..incr as usize).rev() {
            // We have poch1(b, x), but bp is small, so we use backwards
            // recursion to obtain poch1(bp, x)
            let binv = 1.0 / (bp + i as f64);
            poch1 = (poch1 - binv) / (1.0 + x * binv);
        }

        if bp == a {
            result.val = poch1;
            result.err = 2.0 * f64::EPSILON * ((incr.abs() + 1) as f64) * result.val.abs();
            result
        } else {
            // We have poch1(bp, x), but a is < -0.5. We therefore use a
            // reflection formula to obtain poch1(a,x).
            let sinpxx = (PI * x).sin() / x;
            let sinpx2 = (0.5 * PI * x).sin();
            let t1 = sinpxx / (PI * b).tan();
            let t2 = 2.0 * sinpx2 * (sinpx2 / x);
            let trig = t1 - t2;

            result.val = poch1 * (1.0 + x * trig) + trig;
            result.err = ((poch1 * x).abs() + 1.0) * f64::EPSILON * (t1.abs() + t2.abs());
            result.err += 2.0 * f64::EPSILON * ((incr.abs() + 1) as f64) * result.val.abs();
            result
        }
    }
}

// Compute log of the Pochhammer symbol for a > 0 and a + x > 0.
fn lnpoch_pos(a: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let absx = x.abs();

    if absx > 0.1 * a || absx * a.max(2.0).ln() > 0.1 {
        if a < GAMMA_XMAX && a + x < GAMMA_XMAX {
            // If we can do it by calculating the gamma functions
            // directly, then that will be more accurate than
            // doing the subtraction of the logs.
            let g1 = gammainv_e(a);
            let g2 = gammainv_e(a + x);
            result.val = -((g2.val / g1.val).ln());
            result.err = g1.err / (g1.val).abs() + g2.err / (g2.val).abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            // Otherwise, we must do the subtraction
            let lg1 = lngamma_e(a);
            let lg2 = lngamma_e(a + x);
            result.val = lg2.val - lg1.val;
            result.err = lg2.err + lg1.err;
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result
        }
    } else if absx < 0.1 * a && a > 15.0 {
        //  Be careful about the implied subtraction.
        // Note that both a+x and and a must be
        // large here since a is not small
        // and x is not relatively large.
        // So we calculate using Stirling for Log[Gamma(z)].
        //
        //   Log[Gamma(a+x)/Gamma(a)] = x(Log[a]-1) + (x+a-1/2)Log[1+x/a]
        //                              + (1/(1+eps)   - 1) / (12 a)
        //                              - (1/(1+eps)^3 - 1) / (360 a^3)
        //                              + (1/(1+eps)^5 - 1) / (1260 a^5)
        //                              - (1/(1+eps)^7 - 1) / (1680 a^7)
        //                              + ...
        let eps = x / a;
        let den = 1.0 + eps;
        let d3 = den * den * den;
        let d5 = d3 * den * den;
        let d7 = d5 * den * den;
        let c1 = -eps / den;
        let c3 = -eps * (3.0 + eps * (3.0 + eps)) / d3;
        let c5 = -eps * (5.0 + eps * (10.0 + eps * (10.0 + eps * (5.0 + eps)))) / d5;
        let c7 = -eps
            * (7.0 + eps * (21.0 + eps * (35.0 + eps * (35.0 + eps * (21.0 + eps * (7.0 + eps))))))
            / d7;
        let p8 = (1.0 + eps).powi(8);
        let c8 = 1.0 / p8 - 1.0; /* these need not   */
        let c9 = 1.0 / (p8 * (1.0 + eps)) - 1.0; /* be very accurate */
        let a4 = a * a * a * a;
        let a6 = a4 * a * a;
        let ser_1 = c1 + c3 / (30.0 * a * a) + c5 / (105.0 * a4) + c7 / (140.0 * a6);
        let ser_2 = c8 / (99.0 * a6 * a * a) - 691.0 / 360_360.0 * c9 / (a6 * a4);
        let ser = (ser_1 + ser_2) / (12.0 * a);

        let ln_1peps = ln_p1_e(eps);
        let term1 = x * (a / E).ln();
        let term2 = (x + a - 0.5) * ln_1peps.val;

        result.val = term1 + term2 + ser;
        result.err = f64::EPSILON * term1.abs();
        result.err += ((x + a - 0.5) * ln_1peps.err).abs();
        result.err += ln_1peps.val.abs() * f64::EPSILON * (x.abs() + a.abs() + 0.5);
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let poch_rel = pochrel_smallx(a, x);
        let eps = x * poch_rel.val;
        result = ln_p1_e(eps);
        result.err = 2.0 * (x * poch_rel.err / (1.0 + eps)).abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

/// Logarithm of Pochhammer (Apell) symbol:
///     log( (a)_x )
/// where (a)_x := Gamma[a + x]/Gamma[a]
/// and a > 0, a + x > 0.
///
/// # Examples
/// ```
/// assert!((lnpoch_e(1.2,3.5).val - 16.806727225014554213).abs() < 1e-10)
/// ```
pub fn lnpoch_e(a: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if a <= 0.0 || a + x <= 0.0 {
        result.code = SpecFunCode::DomainErr;
        result
    } else if x == 0.0 {
        result
    } else {
        lnpoch_pos(a, x)
    }
}

pub fn lnpoch_sgn_e(a: f64, x: f64) -> (SpecFunResult<f64>, f64) {
    let mut result = SpecFunResult::default();
    if x == 0.0 {
        (result, 1.0)
    } else if a > 0.0 && a + x > 0.0 {
        (lnpoch_e(a, x), 1.0)
    } else if a <= 0.0 && a == a.floor() {
        // Special case for infinite denominator gamma(a)
        if a + x < 0.0 && x == x.floor() {
            let result_pos = lnpoch_pos(-a, -x);
            let f = (a / (a + x)).ln();
            let s = if x % 2.0 == 0.0 { 1.0 } else { -1.0 };
            result.val = f - result_pos.val;
            result.err = result_pos.err + 2.0 * f64::EPSILON * f;
            (result, s)
        } else if a + x == 0.0 {
            // Handle a+x = 0, i.e. gamma(0) / gamma(a)
            // poch(-a,a) == (-1)^a gamma(a+1)
            let (res, sgn) = lngamma_sgn_e(-a + 1.0);
            let s = if -a % 2.0 == 0.0 { 1.0 } else { -1.0 };
            (res, sgn * s)
        } else {
            // Handle finite numberator, Gamma(a+x) for a+x != 0 or neg int
            result.val = f64::NEG_INFINITY;
            result.err = 0.0;
            (result, 1.0)
        }
    } else if a < 0.0 && a + x < 0.0 {
        // Reduce to positive case using reflection.
        let sin_1 = (std::f64::consts::PI * (1.0 - a)).sin();
        let sin_2 = (std::f64::consts::PI * (1.0 - a - x)).sin();

        if sin_1 == 0.0 || sin_2 == 0.0 {
            result.code = SpecFunCode::DomainErr;
            result.val = f64::NAN;
            result.err = f64::NAN;
            (result, 0.0)
        } else {
            let lnp_pos = lnpoch_pos(1.0 - a, -x);
            let lnterm = (sin_1 / sin_2).abs().ln();
            result.val = lnterm - lnp_pos.val;
            result.err = lnp_pos.err;
            result.err +=
                2.0 * f64::EPSILON * ((1.0 - a).abs() + (1.0 - a - a).abs()) * lnterm.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            (result, 1f64.copysign(sin_1 * sin_2))
        }
    } else {
        // Evaluate gamma ratio directly
        let (lg_apn, s_apn) = lngamma_sgn_e(a + x);
        let (lg_a, s_a) = lngamma_sgn_e(a);
        if lg_apn.code == SpecFunCode::Success && lg_a.code == SpecFunCode::Success {
            result.val = lg_apn.val - lg_a.val;
            result.val = lg_apn.val - lg_a.val;
            result.err = lg_apn.err + lg_a.err;
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            (result, s_a * s_apn)
        } else if lg_apn.code == SpecFunCode::DomainErr || lg_a.code == SpecFunCode::DomainErr {
            result.val = f64::NAN;
            result.err = f64::NAN;
            result.code = SpecFunCode::DomainErr;
            (result, 0.0)
        } else {
            result.val = 0.0;
            result.err = 0.0;
            // TODO: Change this to Failure!
            result.code = SpecFunCode::DomainErr;
            (result, 0.0)
        }
    }
}

/// Compute the Pochhammer symbol (x)_a = Gamma(x+a)/Gamma(x).
pub fn poch_e(a: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::default();
    if x == 0.0 {
        result.val = 1.0;
        result
    } else {
        let (lnpoch, sgn) = lnpoch_sgn_e(a, x);
        if lnpoch.val == f64::NEG_INFINITY {
            result
        } else {
            let res = exp_err_e(lnpoch.val, lnpoch.err);
            result.val = res.val * sgn;
            result.err = res.err + 2.0 * f64::EPSILON * result.val.abs();
            result
        }
    }
}

/// Compute the relative Pochhammer symbol [(x)_a - 1] / x.
pub fn pochrel_e(a: f64, x: f64) -> SpecFunResult<f64> {
    let absx = x.abs();
    let absa = a.abs();

    if absx > 0.1 * absa || absx * absa.max(2.0).ln() > 0.1 {
        let mut result = SpecFunResult::default();
        let (lnpoch, sgn) = lnpoch_sgn_e(a, x);
        if lnpoch.val > crate::consts::LN_DBL_MAX {
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.code = SpecFunCode::OverflowErr;
            result
        } else {
            let el = lnpoch.val.exp();
            result.val = (sgn * el - 1.0) / x;
            result.err = result.val.abs() * (lnpoch.err + 2.0 * f64::EPSILON);
            result.err += 2.0 * f64::EPSILON * ((sgn * el).abs() + 1.0) / x.abs();
            result
        }
    } else {
        pochrel_smallx(a, x)
    }
}
