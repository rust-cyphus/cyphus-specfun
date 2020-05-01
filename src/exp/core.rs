use crate::consts::{
    LN_DBL_EPS, LN_DBL_MAX, LN_DBL_MIN, ROOT3_DBL_EPS, SQRT_DBL_MAX, SQRT_DBL_MIN,
};
use crate::gamma::utils::lnfact_int_e;
use crate::result::{SpecFunCode, SpecFunResult, SpecFunResultE10};
use std::f64::consts::{LN_10};

pub(crate) fn exp_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    if x > LN_DBL_MAX {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result
    } else if x < LN_DBL_MIN {
        result.code = SpecFunCode::UnderflowErr;
        result
    } else {
        result.val = x.exp();
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

/// Compute exp(x)
pub(crate) fn exp_e10_e(x: f64) -> SpecFunResultE10<f64> {
    let mut result = SpecFunResultE10 {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
        e10: 0,
    };

    if x > (i32::MAX - 1) as f64 {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("exp_e10_e", &[x]);
        result
    } else if x < (i32::MIN + 1) as f64 {
        result.val = 0.0;
        result.err = 0.0;
        result.code = SpecFunCode::UnderflowErr;
        result.issue_warning("exp_e10_e", &[x]);
        result
    } else {
        let n: i32 = if x > LN_DBL_MAX || x < LN_DBL_MIN {
            (x / LN_10).floor() as i32
        } else {
            0
        };
        result.val = (x - n as f64 * LN_10).exp();
        result.err = 2.0 * (1.0 + x.abs()) * f64::EPSILON * result.val.abs();
        result.e10 = n;
        result
    }
}

pub(crate) fn exp_mult_e(x: f64, y: f64) -> SpecFunResult<f64> {
    let ay = y.abs();

    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if ay < f64::EPSILON {
        result
    } else if (x < 0.5 * LN_DBL_MAX && x > 0.5 * LN_DBL_MIN)
        && ay < 0.8 * SQRT_DBL_MAX
        && ay > 1.2 * SQRT_DBL_MIN
    {
        let ex = x.exp();
        result.val = y * ex;
        result.err = (2.0 + x.abs()) * f64::EPSILON * result.val.abs();
        result
    } else {
        let ly = ay.ln();
        let lnr = x + ly;

        if lnr > LN_DBL_MAX - 0.01 {
            result.code = SpecFunCode::OverflowErr;
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else if lnr < LN_DBL_MIN + 0.01 {
            result.code = SpecFunCode::UnderflowErr;
            result.val = 0.0;
            result.err = 0.0;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else {
            let sy = y.signum();
            let m = x.floor();
            let n = ly.floor();
            let a = x - m;
            let b = ly - n;
            let berr = 2.0 * f64::EPSILON * (ly.abs() + n.abs());

            result.val = sy * (m + n).exp() * (a + b).exp();
            result.err = berr * result.val.abs();
            result.err += 2.0 * f64::EPSILON * (m + n + 1.0) * result.val.abs();
            result
        }
    }
}

pub(crate) fn exp_mult_e10_e(x: f64, y: f64) -> SpecFunResultE10<f64> {
    let ay = y.abs();

    let mut result = SpecFunResultE10 {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
        e10: 0,
    };

    if ay < f64::EPSILON {
        result
    } else if (x < 0.5 * LN_DBL_MAX && x > 0.5 * LN_DBL_MIN)
        && ay < 0.8 * SQRT_DBL_MAX
        && ay > 1.2 * SQRT_DBL_MIN
    {
        let ex = x.exp();
        result.val = y * ex;
        result.err = (2.0 + x.abs()) * f64::EPSILON * result.val.abs();
        result
    } else {
        let ly = ay.abs();
        let l10_val = (x + ly) / LN_10;

        if l10_val > (i32::MAX - 1) as f64 {
            result.code = SpecFunCode::OverflowErr;
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else if l10_val < (i32::MIN + 1) as f64 {
            result.code = SpecFunCode::UnderflowErr;
            result.val = 0.0;
            result.err = 0.0;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else {
            let sy = y.signum();
            let n = l10_val.floor() as i32;
            let arg_val = (l10_val - n as f64) * LN_10;
            let arg_err = 2.0 * f64::EPSILON * (x.abs() + ly.abs() + LN_10 * (n as f64).abs());

            result.val = sy * arg_val.exp();
            result.err = arg_err * result.val.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result.e10 = n;
            result
        }
    }
}

/// Exponentiate and multiply by a given factor: y * Exp(x)
pub(crate) fn exp_mult_err_e(x: f64, dx: f64, y: f64, dy: f64) -> SpecFunResult<f64> {
    let ay = y.abs();

    if ay < f64::EPSILON {
        SpecFunResult {
            val: 0.0,
            err: (dy * x.exp()).abs(),
            code: SpecFunCode::Success,
        }
    } else if (x < 0.5 * f64::MAX.ln() && x > 0.5 * f64::MIN_POSITIVE.ln())
        && (ay < 0.8 * f64::MAX.sqrt() && ay > 1.0 * f64::MIN_POSITIVE.sqrt())
    {
        let ex = x.exp();
        let val = y * ex;
        let mut err = ex * (dy.abs() + (y * dx).abs());
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let ly = ay.ln();
        let lnr = x + ly;

        if lnr > f64::MAX.ln() - 0.01 {
            let result = SpecFunResult {
                val: f64::INFINITY,
                err: f64::INFINITY,
                code: SpecFunCode::OverflowErr,
            };
            result.issue_warning("exp_mult_err_e", &[x, dx, y, dy]);
            result
        } else if lnr < f64::MIN_POSITIVE.ln() + 0.01 {
            let result = SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::UnderflowErr,
            };
            result.issue_warning("exp_mult_err_e", &[x, dx, y, dy]);
            result
        } else {
            let sy = y.signum();
            let M = x.floor();
            let N = ly.floor();
            let a = x - M;
            let b = ly - N;
            let eMN = (M + N).exp();
            let eab = (a + b).exp();
            let val = sy * eMN * eab;
            let mut err = eMN * eab * 2.0 * f64::EPSILON;
            err += eMN * eab * (dy / y).abs();
            err += eMN * eab * dx.abs();
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        }
    }
}

pub(crate) fn exp_mult_err_e10_e(x: f64, dx: f64, y: f64, dy: f64) -> SpecFunResultE10<f64> {
    let mut result = SpecFunResultE10 {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
        e10: 0,
    };

    let ay = y.abs();

    if ay < f64::EPSILON {
        result.err = (dy * x.exp()).abs();
        result
    } else if (x < 0.5 * LN_DBL_MAX && x > 0.5 * LN_DBL_MIN)
        && ay < 0.8 * SQRT_DBL_MAX
        && ay > 1.2 * SQRT_DBL_MIN
    {
        let ex = x.exp();
        result.val = y * ex;
        result.err = ex * (dy.abs() + (y * dx).abs());
        result.err += (2.0 + x.abs()) * f64::EPSILON * result.val.abs();
        result
    } else {
        let ly = ay.abs();
        let l10_val = (x + ly) / LN_10;

        if l10_val > (i32::MAX - 1) as f64 {
            result.code = SpecFunCode::OverflowErr;
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else if l10_val < (i32::MIN + 1) as f64 {
            result.code = SpecFunCode::UnderflowErr;
            result.val = 0.0;
            result.err = 0.0;
            result.issue_warning("exp_mult_err_e10_e", &[x]);
            result
        } else {
            let sy = y.signum();
            let n = l10_val.floor() as i32;
            let arg_val = (l10_val - n as f64) * LN_10;
            let arg_err = dy / y.abs() + dx + 2.0 * f64::EPSILON * arg_val.abs();

            result.val = sy * arg_val.exp();
            result.err = arg_err * result.val.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result.e10 = n;
            result
        }
    }
}

pub(crate) fn expm1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let cut = 0.002;

    if x < LN_DBL_MIN {
        result.val = -1.0;
        result.err = f64::EPSILON;
    } else if x < -cut {
        result.val = x.exp() - 1.0;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < cut {
        result.val = x * 0.2_f64
            .mul_add(x, 1.0)
            .mul_add(0.25 * x, 1.0)
            .mul_add(x / 3.0, 1.0)
            .mul_add(0.5 * x, 1.0);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < LN_DBL_MAX {
        result.val = x.exp() - 1.0;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("expm1_e", &[x]);
    }

    result
}

pub(crate) fn exprel_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let cut = 0.002;

    if x < LN_DBL_MIN {
        result.val = -1.0 / x;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < -cut {
        result.val = (x.exp() - 1.0) / x;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < cut {
        result.val = 0.2_f64
            .mul_add(x, 1.0)
            .mul_add(0.25 * x, 1.0)
            .mul_add(x / 3.0, 1.0)
            .mul_add(0.5 * x, 1.0);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < LN_DBL_MAX {
        result.val = (x.exp() - 1.0) / x;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("exprel_e", &[x]);
    }

    result
}

pub(crate) fn exprel_2_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let cut = 0.002;

    if x < LN_DBL_MIN {
        result.val = -2.0 / x * (1.0 + 1.0 / x);
        result.err = f64::EPSILON * result.val.abs();
    } else if x < -cut {
        result.val = 2.0 * (x.exp() - 1.0 - x) / (x * x);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < cut {
        result.val = 6.0_f64
            .recip()
            .mul_add(x, 1.0)
            .mul_add(0.2 * x, 1.0)
            .mul_add(0.25 * x, 1.0)
            .mul_add(x / 3.0, 1.0);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else if x < LN_DBL_MAX {
        result.val = 2.0 * (x.exp() - 1.0 - x) / (x * x);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("exprel_2_e", &[x]);
    }

    result
}

/// Evaluate continued fraction for exprel
/// Ref: Abramowitz + Stegun, 4.2.41
pub(crate) fn exprel_n_CF_e(nn: usize, x: f64) -> SpecFunResult<f64> {
    let RECUR_BIG = SQRT_DBL_MAX;
    let maxiter = 5000;
    let mut n = 1;
    let mut Anm2: f64 = 1.0;
    let mut Bnm2: f64 = 0.0;
    let mut Anm1: f64 = 0.0;
    let mut Bnm1: f64 = 1.0;
    let a1: f64 = 1.0;
    let b1: f64 = 1.0;
    let a2: f64 = -x;
    let b2: f64 = (nn + 1) as f64;

    let mut An = b1 * Anm1 + a1 * Anm2; /* A1 */
    let mut Bn = b1 * Bnm1 + a1 * Bnm2; /* B1 */

    // One explicit step, before we get to the main pattern.
    n += 1;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    An = b2 * Anm1 + a2 * Anm2; /* A2 */
    Bn = b2 * Bnm1 + a2 * Bnm2; /* B2 */

    let mut ffn = An / Bn;

    while n < maxiter {
        n += 1;

        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        let an = if n % 2 == 1 {
            ((n - 1) / 2) as f64 * x
        } else {
            -((nn + (n / 2) - 1) as f64) * x
        };
        let bn = (nn + n - 1) as f64;
        An = bn * Anm1 + an * Anm2;
        Bn = bn * Bnm1 + an * Bnm2;

        if An.abs() > RECUR_BIG || Bn.abs() > RECUR_BIG {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
            Bnm2 /= RECUR_BIG;
        }

        let old_fn = ffn;
        ffn = An / Bn;
        let del = old_fn / ffn;

        if (del - 1.0).abs() < 2.0 * f64::EPSILON {
            break;
        }
    }

    let mut result = SpecFunResult {
        val: ffn,
        err: 4.0 * (n + 1) as f64 * f64::EPSILON * ffn.abs(),
        code: SpecFunCode::Success,
    };
    if n == maxiter {
        result.code = SpecFunCode::MaxIterErr;
        result.issue_warning("", &[nn as f64, x]);
    }

    result
}

pub(crate) fn exprel_n_e(n: i32, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    let nd = n as f64;

    if n < 0 {
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("exprel_n_e", &[n as f64, x]);
        result
    } else if x.abs() <= f64::EPSILON {
        result.val = 1.0;
        result.err = 0.0;
        result
    } else if x.abs() < ROOT3_DBL_EPS * n as f64 {
        result.val = 1.0 + x / ((n + 1) as f64) * (1.0 + x / ((n + 2) as f64));
        result.err = 2.0 * f64::EPSILON;
        result
    } else if n == 0 {
        exp_e(x)
    } else if n == 1 {
        exprel_e(x)
    } else if n == 2 {
        exprel_2_e(x)
    } else if x > nd && (-x + nd * 1.0 + (x / nd).ln()) < LN_DBL_EPS {
        // x is much larger than n.
        // Ignore polynomial part, so
        // exprel_N(x) ~= e^x N!/x^N
        let lnfn = lnfact_int_e(n as usize);
        let lnterm = nd * x.ln();
        let lnr_val = x + lnfn.val - lnterm;
        let lnr_err = (f64::EPSILON * (x.abs() + lnfn.val.abs() + lnterm.abs())) + lnfn.err;
        exp_err_e(lnr_val, lnr_err)
    } else if x > nd {
        // Write the identity
        // exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
        // then use the asymptotic expansion
        // Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)

        let lnx = x.ln();
        let lnfn = lnfact_int_e(n as usize);
        let lgn = lnfn.val - nd.ln();
        let lnpre_val = x + lnfn.val - nd * lnx;
        let lnpre_err = f64::EPSILON * (x.abs() + lnfn.val.abs() + (nd * lnx).abs()) + lnfn.err;
        if lnpre_val < LN_DBL_MAX - 5.0 {
            let pre = exp_err_e(lnpre_val, lnpre_err);
            let ln_bigg_ratio_pre = -x + (nd - 1.0) * lnx - lgn;
            let mut biggsum = 1.0;
            let mut term = 1.0;
            for k in 1..n {
                term *= (n - k) as f64 / x;
                biggsum += term;
            }
            let bigg_ratio = exp_mult_e(ln_bigg_ratio_pre, biggsum);
            if bigg_ratio.code == SpecFunCode::Success {
                result.val = pre.val * (1.0 - bigg_ratio.val);
                result.err = pre.val * (2.0 * f64::EPSILON + bigg_ratio.err);
                result.err += pre.err * (1.0 - bigg_ratio.val).abs();
                result.err += 2.0 * f64::EPSILON * result.val.abs();
                result
            } else {
                result.val = 0.0;
                result.err = 0.0;
                result
            }
        } else {
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.code = SpecFunCode::OverflowErr;
            result
        }
    } else if x > -10.0 * nd {
        exprel_n_CF_e(n as usize, x)
    } else {
        // x -> -Inf asymptotic:
        // exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
        //             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)

        let mut sum = 1.0;
        let mut term = 1.0;
        for k in 1..n {
            term *= (n - k) as f64 / x;
            sum += term;
        }
        result.val = -nd / x * sum;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

/// Exponentiate a quantity with an associated error.
pub(crate) fn exp_err_e(x: f64, dx: f64) -> SpecFunResult<f64> {
    let adx = dx.abs();

    if x + adx > LN_DBL_MAX {
        let result = SpecFunResult {
            val: f64::INFINITY,
            err: f64::INFINITY,
            code: SpecFunCode::OverflowErr,
        };
        result.issue_warning("exp_err_e", &[x, dx]);
        result
    } else if x - adx < LN_DBL_MIN {
        let result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::UnderflowErr,
        };
        result.issue_warning("exp_err_e", &[x, dx]);
        result
    } else {
        let ex = x.exp();
        let edx = adx.exp();
        let val = ex;
        let mut err = ex * f64::EPSILON.min(edx - 1.0 / edx);
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

pub(crate) fn exp_err_e10_e(x: f64, dx: f64) -> SpecFunResultE10<f64> {
    let adx = dx.abs();
    let mut result = SpecFunResultE10 {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
        e10: 0,
    };

    if x + adx > (i32::MAX - 1) as f64 {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("exp_err_e10_e", &[x, dx]);
    } else if x - adx < (i32::MIN + 1) as f64 {
        result.val = 0.0;
        result.err = 0.0;
        result.code = SpecFunCode::UnderflowErr;
        result.issue_warning("exp_err_e10_e", &[x, dx]);
    } else {
        let n = (x / LN_10).floor();
        let ex = (x - n * LN_10).exp();
        result.val = ex;
        result.err = ex * (2.0 * f64::EPSILON * (x.abs() + 1.0) + adx);
        result.e10 = n as i32;
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::consts::{LN_DBL_MAX, SQRT_DLB_EPS};
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_exp_e() {
        test_sf_check_result(exp_e(-10.0), (-10.0_f64).exp(), TOL0);
        test_sf_check_result(exp_e(10.0), (10.0_f64).exp(), TOL0);
    }
    #[test]
    fn test_exp_e10_e() {
        test_sf_e10(exp_e10_e(1.0), std::f64::consts::E, 0, TOL0);
        assert!(exp_e10_e(1.0).err <= TOL1, "Error too large.");
        test_sf_e10(exp_e10_e(2000.0), 3.881_180_194_283_637_3, 868, TOL3);
        assert!(exp_e10_e(2000.0).err <= TOL5, "Error too large.");
    }
    #[test]
    fn test_exp_err_e() {
        test_sf_check_result(exp_err_e(-10.0, TOL1), (-10.0_f64).exp(), TOL1);
        test_sf_check_result(exp_err_e(10.0, TOL1), (10.0_f64).exp(), TOL1);
    }
    #[test]
    fn test_exp_err_e10_e() {
        test_sf_e10(exp_err_e10_e(1.0, SQRT_TOL0), std::f64::consts::E, 0, TOL1);
        assert!(exp_e10_e(1.0).err <= 32.0 * SQRT_TOL0, "Error too large.");

        test_sf_e10(exp_err_e10_e(2000.0, 1e-10), 3.881_180_194_283_637_3, 868, TOL3);
        assert!(exp_e10_e(2000.0).err <= 1e-7, "Error too large.");
    }
    #[test]
    fn test_exp_mult_e() {
        test_sf_check_result_and_code(
            exp_mult_e(-10.0, 1e-6),
            1e-6 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(-10.0, 2.0),
            2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(-10.0, -2.0),
            -2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(10.0, 1e-6),
            1e-6 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(10.0, 2.0),
            2.0 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success,
        );

        let x = 0.8 * LN_DBL_MAX;
        test_sf_check_result_and_code(
            exp_mult_e(x, 1.00001),
            1.00001 * x.exp(),
            TOL3,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(x, 1.000001),
            1.000001 * x.exp(),
            TOL3,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(x, 100.0),
            100.0 * x.exp(),
            TOL3,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(x, 1e20),
            1e20 * x.exp(),
            TOL3,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            exp_mult_e(x, (-x).exp() * LN_2.exp()),
            2.0,
            TOL4,
            SpecFunCode::Success,
        );
    }
}
