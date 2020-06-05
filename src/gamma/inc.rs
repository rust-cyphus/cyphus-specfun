use super::mono::{gammastar_e, lngamma_e};
use crate::exp::core::exprel_n_cf_e;
use crate::logarithm::ln_p1_mx_e;
use crate::result::{SpecFunCode, SpecFunResult};

// The dominant part,
// D(a,x) := x^a e^(-x) / Gamma(a+1)
#[allow(dead_code)]
fn gamma_inc_d(a: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if a < 10.0 {
        let lg = lngamma_e(a + 1.0);
        let lnr = a * x.ln() - x - lg.val;
        result.val = lnr.exp();
        result.err = 2.0 * f64::EPSILON * (lnr.abs() + 1.0) * result.val.abs();
    } else {
        let mut ln_term = SpecFunResult::<f64>::default();
        if x < 0.5 * a {
            let u = x / a;
            let ln_u = u.ln();
            ln_term.val = ln_u - u + 1.0;
            ln_term.err = (ln_u.abs() + u.abs() + 1.0) * f64::EPSILON;
        } else {
            let mu = (x - a) / a;
            ln_term = ln_p1_mx_e(mu); /* log(1+mu) - mu */
            /* Propagate cancellation error from x-a, since the absolute
            error of mu=x-a is DBL_EPSILON */
            ln_term.err += f64::EPSILON * mu.abs();
        };
        let gstar = gammastar_e(a);
        let term1 = (a * ln_term.val).exp() / (2.0 * std::f64::consts::PI * a).sqrt();
        result.val = term1 / gstar.val;
        result.err = 2.0 * f64::EPSILON * ((a * ln_term.val).abs() + 1.0) * result.val.abs();
        /* Include propagated error from log term */
        result.err += a.abs() * ln_term.err * result.val.abs();
        result.err += gstar.err / gstar.val.abs() * result.val.abs();
    }
    result
}

// P series representation.
#[allow(dead_code)]
fn gamma_inc_p_series(a: f64, x: f64) -> SpecFunResult<f64> {
    let nmax = 10000;
    let mut result = SpecFunResult::<f64>::default();

    let d = gamma_inc_d(a, x);

    // Approximating the terms of the series using Stirling's
    // approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
    // convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
    // -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
    // below detects cases where the minimum value of n is > 5000

    if x > 0.995 * a && a > 1e5 {
        // Difficult case: try continued fraction
        let cf_res = exprel_n_cf_e(a, x);
        result.val = d.val * cf_res.val;
        result.err = (d.val * cf_res.err).abs() + (d.err * cf_res.val).abs();
        return result;
    }

    // Series would require excessive number of terms
    if x > a + nmax as f64 {
        result.code = SpecFunCode::MaxIterErr;
        result.issue_warning("gamma_inc_p_series_e", &[a, x]);
        return result;
    }

    // Normal case: sum the series
    let mut sum = 1.0;
    let mut term = 1.0;

    // Handle lower part of the series where t_n is increasing, |x| > a+n

    let nlow = if x > a { (x - a).ceil() as usize } else { 0 };

    for n in 1..nlow {
        term *= x / (a + n as f64);
        sum += term;
    }

    // Handle upper part of the series where t_n is decreasing, |x| < a+n
    let mut n = nlow - 1;
    while n < nmax {
        term *= x / (a + n as f64);
        sum += term;
        if (term / sum).abs() < f64::EPSILON {
            break;
        }
        n += 1;
    }

    //  Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1))
    let remainder = {
        let tnp1 = (x / (a + n as f64)) * term;
        tnp1 / (1.0 - x / (a + n as f64 + 1.0))
    };

    result.val = d.val * sum;
    result.err = d.err * sum.abs() + (d.val * remainder).abs();
    result.err += (1.0 + n as f64) * f64::EPSILON * result.val.abs();

    if n == nmax && (remainder / sum).abs() > crate::consts::SQRT_DLB_EPS {
        result.code = SpecFunCode::MaxIterErr;
        result.issue_warning("gamma_inc_p_series_e", &[a, x]);
    }

    result
}

// Q large x asymptotic
#[allow(dead_code)]
fn gamma_inc_q_large_x(a: f64, x: f64) -> SpecFunResult<f64> {
    let nmax = 5000;
    let mut result = SpecFunResult::<f64>::default();

    let d = gamma_inc_d(a, x);

    let mut sum = 1.0;
    let mut term = 1.0;
    let mut last = 1.0;

    for n in 1..nmax {
        term *= (a - n as f64) / x;
        if (term / last).abs() > 1.0 {
            break;
        }
        if (term / sum).abs() < f64::EPSILON {
            break;
        }
        sum += term;
        last = term;

        if n == nmax - 1 {
            result.code = SpecFunCode::MaxIterErr;
            result.issue_warning("gamma_inc_q_large_x", &[a, x]);
        }
    }

    result.val = d.val * (a / x) * sum;
    result.err = d.err * ((a / x) * sum).abs();
    result.err += 2.0 * f64::EPSILON * (result.val).abs();

    result
}
