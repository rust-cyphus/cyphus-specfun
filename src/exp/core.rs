use crate::consts::{
    LN_DBL_EPS, LN_DBL_MAX, LN_DBL_MIN, ROOT3_DBL_EPS, SQRT_DBL_MAX, SQRT_DBL_MIN,
};
use crate::gamma::mono::lnfact_e;
use crate::result::{SpecFunCode, SpecFunResult, SpecFunResultE10};
use std::f64::consts::LN_10;

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

    if y == 0.0 {
        result
    } else if (x < 0.5 * LN_DBL_MAX && x > 0.5 * LN_DBL_MIN)
        && (ay < 0.8 * SQRT_DBL_MAX && ay > 1.2 * SQRT_DBL_MIN)
    {
        result.val = y * x.exp();
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

    if y == 0.0 {
        result
    } else if (x < 0.5 * LN_DBL_MAX && x > 0.5 * LN_DBL_MIN)
        && (ay < 0.8 * SQRT_DBL_MAX && ay > 1.2 * SQRT_DBL_MIN)
    {
        result.val = y * x.exp();
        result.err = (2.0 + x.abs()) * f64::EPSILON * result.val.abs();
        result
    } else {
        let ly = ay.ln();
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
            let n = l10_val.floor();
            let arg_val = (l10_val - n) * LN_10;
            let arg_err = 2.0 * f64::EPSILON * (x.abs() + ly.abs() + LN_10 * n.abs());

            result.val = sy * arg_val.exp();
            result.err = arg_err * result.val.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result.e10 = n as i32;
            result
        }
    }
}

/// Exponentiate and multiply by a given factor: y * Exp(x)
pub(crate) fn exp_mult_err_e(x: f64, dx: f64, y: f64, dy: f64) -> SpecFunResult<f64> {
    let ay = y.abs();

    if y == 0.0 {
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
            let m = x.floor();
            let n = ly.floor();
            let a = x - m;
            let b = ly - n;
            let emn = (m + n).exp();
            let eab = (a + b).exp();
            let val = sy * emn * eab;
            let mut err = emn * eab * 2.0 * f64::EPSILON;
            err += emn * eab * (dy / y).abs();
            err += emn * eab * dx.abs();
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

    if y == 0.0 {
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
        let ly = ay.ln();
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
pub(crate) fn exprel_n_cf_e(nn: usize, x: f64) -> SpecFunResult<f64> {
    let recur_big = SQRT_DBL_MAX;
    let maxiter = 5000;
    let mut n = 1;
    let mut anm2: f64 = 1.0;
    let mut bnm2: f64 = 0.0;
    let mut anm1: f64 = 0.0;
    let mut bnm1: f64 = 1.0;
    let a1: f64 = 1.0;
    let b1: f64 = 1.0;
    let a2: f64 = -x;
    let b2: f64 = (nn + 1) as f64;

    let mut aan = b1 * anm1 + a1 * anm2; /* A1 */
    let mut bbn = b1 * bnm1 + a1 * bnm2; /* B1 */

    // One explicit step, before we get to the main pattern.
    n += 1;
    anm2 = anm1;
    bnm2 = bnm1;
    anm1 = aan;
    bnm1 = bbn;
    aan = b2 * anm1 + a2 * anm2; /* A2 */
    bbn = b2 * bnm1 + a2 * bnm2; /* B2 */

    let mut ffn = aan / bbn;

    while n < maxiter {
        n += 1;

        anm2 = anm1;
        bnm2 = bnm1;
        anm1 = aan;
        bnm1 = bbn;
        let an = if n % 2 == 1 {
            ((n - 1) / 2) as f64 * x
        } else {
            -((nn + (n / 2) - 1) as f64) * x
        };
        let bn = (nn + n - 1) as f64;
        aan = bn * anm1 + an * anm2;
        bbn = bn * bnm1 + an * bnm2;

        if aan.abs() > recur_big || bbn.abs() > recur_big {
            aan /= recur_big;
            bbn /= recur_big;
            anm1 /= recur_big;
            bnm1 /= recur_big;
        }

        let old_fn = ffn;
        ffn = aan / bbn;
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

pub(crate) fn exprel_n_e(n: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    let nd = n as f64;

    if x.abs() <= f64::MIN_POSITIVE {
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
    } else if x > nd && -x + nd * (1.0 + (x / nd).ln()) < LN_DBL_EPS {
        // x is much larger than n.
        // Ignore polynomial part, so
        // exprel_N(x) ~= e^x N!/x^N
        let lnfn = lnfact_e(n as usize);
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
        let lnfn = lnfact_e(n);
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
        exprel_n_cf_e(n, x)
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
