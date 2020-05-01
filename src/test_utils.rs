use crate::result::{SpecFunCode, SpecFunResult, SpecFunResultE10};

pub(crate) enum TestResult {
    InCons = 1,
    ErrNeg = 2,
    TolBad = 4,
    RetBad = 8,
    ErrBad = 16,
    ErrBig = 32,
    ExpBad = 64,
}

pub(crate) fn test_sf_frac_diff(x1: f64, x2: f64) -> f64 {
    if x1.abs() < f64::EPSILON && x2.abs() < f64::EPSILON {
        0.0
    } else if x1.abs() < f64::EPSILON {
        x2.abs()
    } else if (x1 <= f64::MAX && x2 <= f64::MAX) && (x1 + x2).abs() > f64::EPSILON {
        ((x1 - x2) / (x1 + x2)).abs()
    } else {
        1.0
    }
}

// Check a result against a given value and tolerance.
pub(crate) fn test_sf_check_result(r: SpecFunResult<f64>, val: f64, tol: f64) {
    let mut s: usize = 0;
    let mut f: f64 = 0.0;
    let mut d: f64;

    if r.val.is_nan() || val.is_nan() {
        s = if r.val.is_nan() != val.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
    } else if r.val.is_infinite() || val.is_infinite() {
        s = if r.val.is_infinite() != val.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
    } else {
        f = test_sf_frac_diff(val, r.val);
        d = (val - r.val).abs();

        if d > 2.0 * r.err {
            s |= TestResult::InCons as usize;
        }
        if r.err < 0.0 {
            s |= TestResult::ErrNeg as usize;
        }
        if r.err.is_infinite() {
            s |= TestResult::ErrBad as usize;
        }

        if f > tol {
            s |= TestResult::TolBad as usize;
        }
    }

    if s != 0 {
        println!("expected: {:?}", val);
        println!(
            "obtained: {:?} +/- {:?} (rel={:?})",
            r.val,
            r.err,
            r.err / (r.val + r.err).abs()
        );
        println!("fracdiff: {:?}", f);
        println!("tolerance: {:?}", tol);
    }

    assert!(
        s & (TestResult::InCons as usize) == 0,
        "value/expected not consistent within reported error"
    );
    assert!(
        s & (TestResult::ErrNeg as usize) == 0,
        "reported error negative"
    );
    assert!(
        s & (TestResult::ErrBad as usize) == 0,
        "reported error is bad"
    );
    assert!(
        s & (TestResult::ErrBig as usize) == 0,
        "reported error is much too big"
    );
    assert!(
        s & (TestResult::TolBad as usize) == 0,
        "value not within tolerance of expected value"
    );
}

pub(crate) fn test_sf_check_result_and_code(
    r: SpecFunResult<f64>,
    val: f64,
    tol: f64,
    code: SpecFunCode,
) {
    let mut s: usize = 0;
    let mut f: f64 = 0.0;
    let mut d: f64;

    if r.val.is_nan() || val.is_nan() {
        s = if r.val.is_nan() != val.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
    } else if r.val.is_infinite() || val.is_infinite() {
        s = if r.val.is_infinite() != val.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
    } else {
        f = test_sf_frac_diff(val, r.val);
        d = (val - r.val).abs();

        if d > 2.0 * r.err {
            s |= TestResult::InCons as usize;
        }
        if r.err < 0.0 {
            s |= TestResult::ErrNeg as usize;
        }
        if r.err.is_infinite() {
            s |= TestResult::ErrBad as usize;
        }

        if f > tol {
            s |= TestResult::TolBad as usize;
        }
    }

    if s != 0 {
        println!("expected: {:?}", val);
        println!(
            "obtained: {:?} +/- {:?} (rel={:?})",
            r.val,
            r.err,
            r.err / (r.val + r.err).abs()
        );
        println!("fracdiff: {:?}", f);
        println!("tolerance: {:?}", tol);
    }

    assert!(
        s & (TestResult::InCons as usize) == 0,
        "value/expected not consistent within reported error"
    );
    assert!(
        s & (TestResult::ErrNeg as usize) == 0,
        "reported error negative"
    );
    assert!(
        s & (TestResult::ErrBad as usize) == 0,
        "reported error is bad"
    );
    assert!(
        s & (TestResult::ErrBig as usize) == 0,
        "reported error is much too big"
    );
    assert!(
        s & (TestResult::TolBad as usize) == 0,
        "value not within tolerance of expected value"
    );
    assert!(r.code == code, "Retcodes are not equal.");
}

pub(crate) fn test_sf_check_e10(e10: i32, e10_in: i32) {
    let mut s = 0;

    if e10 != e10_in {
        s = TestResult::ExpBad as usize;
    }

    if s != 0 {
        println!("Expected exponent: 10^{:?}", e10_in);
        println!("Obtained exponent: 10^{:?}", e10);
    }

    assert!(
        s & TestResult::ExpBad as usize == 0,
        "Exponent is incorrect."
    );
}

pub(crate) fn test_sf_check_val(rval: f64, val: f64, tol: f64) {
    let mut s: usize = 0;
    let f = test_sf_frac_diff(val, rval);

    if f > tol {
        s |= TestResult::TolBad as usize;
    }

    if s != 0 {
        println!("expected: {:?}", val);
        println!("obtained:{:?}", rval);
        println!("fracdiff: {:?}", f);
    }
    assert!(
        s & TestResult::TolBad as usize == 0,
        "value not within tolerance of expected value"
    );
}

pub(crate) fn test_sf_check_result_relax(r: SpecFunResult<f64>, val: f64, tol: f64) {
    let mut s = 0;
    let f = test_sf_frac_diff(val, r.val);

    if f > tol.max(1.0) {
        s |= TestResult::InCons as usize;
    }
    if r.err < 0.0 {
        s |= TestResult::ErrNeg as usize;
    }
    if r.err.is_infinite() {
        s |= TestResult::ErrBad as usize;
    }
    if f > tol {
        s |= TestResult::ErrBad as usize;
    }

    if s != 0 {
        println!("expected: {:?}", val);
        println!(
            "obtained: {:?} +/- {:?} (rel={:?})",
            r.val,
            r.err,
            r.err / (r.val + r.err).abs()
        );
        println!("fracdiff: {:?}", f);
    }
    assert!(
        s & (TestResult::InCons as usize) == 0,
        "value/expected not consistent within reported error"
    );
    assert!(
        s & (TestResult::ErrNeg as usize) == 0,
        "reported error negative"
    );
    assert!(
        s & (TestResult::ErrBad as usize) == 0,
        "reported error is bad"
    );
    assert!(
        s & (TestResult::TolBad as usize) == 0,
        "value not within tolerance of expected value"
    );
}

pub(crate) fn test_sf_check_return(val_return: usize, val_expected: usize) {
    assert!(
        val_return == val_expected,
        format!("unexpected return code: {:?}", val_return)
    );
}

pub(crate) fn test_sf_e10(re: SpecFunResultE10<f64>, val_in: f64, e10_in: i32, tol: f64) {
    let r = SpecFunResult {
        val: re.val,
        err: re.err,
        code: re.code,
    };
    test_sf_check_result(r, val_in, tol);
    test_sf_check_e10(re.e10, e10_in);
}
