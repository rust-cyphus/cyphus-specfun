#[allow(dead_code)]
pub(crate) enum TestResult {
    InCons = 1,
    ErrNeg = 2,
    TolBad = 4,
    RetBad = 8,
    ErrBad = 16,
    ErrBig = 32,
    ExpBad = 64,
}

#[allow(dead_code)]
pub(crate) fn test_sf_frac_diff(x1: f64, x2: f64) -> f64 {
    if x1 == 0.0 && x2 == 0.0 {
        0.0
    } else if x1 == 0.0 {
        x2.abs()
    } else if (x1 <= f64::MAX && x2 <= f64::MAX) && x1 + x2 != 0.0 {
        ((x1 - x2) / (x1 + x2)).abs()
    } else {
        1.0
    }
}

#[macro_export]
macro_rules! test_check_result_and_code {
    ($func:ident, ($($arg:expr ),*), $val:expr, $tol:expr, $code:expr) => {
        let r = $func($($arg),*);
        let mut s:usize = 0;
        let mut f: f64 = 0.0;
        let mut d: f64 = 0.0;
        let ival :f64 = ($val as f64);

        if r.val.is_nan() ||ival.is_nan() {
            s = if r.val.is_nan() != ival.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else if r.val.is_infinite() || ival.is_infinite() {
            s = if r.val.is_infinite() != ival.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else {
            f = test_sf_frac_diff(ival, r.val);
            d = (ival - r.val).abs();

            if d > 2.5 * r.err {
                s |= TestResult::InCons as usize;
            }
            if r.err < 0.0 {
                s |= TestResult::ErrNeg as usize;
            }
            if r.err.is_infinite() {
                s |= TestResult::ErrBad as usize;
            }

            if f > $tol {
                s |= TestResult::TolBad as usize;
            }
        }

        if s != 0 {
            println!("expected: {:?}", ival);
            println!(
                "obtained: {:?} +/- {:?} (rel={:?})",
                r.val,
                r.err,
                r.err / (r.val + r.err).abs()
            );
            println!("fracdiff: {:?}", f);
            println!("tolerance: {:?}", $tol);
            println!("absdiff: {:?}", d);
            println!("call: {:?}({:?})",stringify!($func),($($arg),*));
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
        assert!(r.code == $code, "Retcodes are not equal.");
    };

    ($res:expr, $val:expr, $tol:expr, $code:expr) =>{
        let r = $res;
        let mut s:usize = 0;
        let mut f: f64 = 0.0;
        let mut d: f64 = 0.0;
        let ival :f64 = ($val as f64);

        if r.val.is_nan() ||ival.is_nan() {
            s = if r.val.is_nan() != ival.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else if r.val.is_infinite() || ival.is_infinite() {
            s = if r.val.is_infinite() != ival.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else {
            f = test_sf_frac_diff(ival, r.val);
            d = (ival - r.val).abs();

            if d > 2.0 * r.err {
                s |= TestResult::InCons as usize;
            }
            if r.err < 0.0 {
                s |= TestResult::ErrNeg as usize;
            }
            if r.err.is_infinite() {
                s |= TestResult::ErrBad as usize;
            }

            if f > $tol {
                s |= TestResult::TolBad as usize;
            }
        }

        if s != 0 {
            println!("expected: {:?}", ival);
            println!(
                "obtained: {:?} +/- {:?} (rel={:?})",
                r.val,
                r.err,
                r.err / (r.val + r.err).abs()
            );
            println!("fracdiff: {:?}", f);
            println!("tolerance: {:?}", $tol);
            println!("absdiff: {:?}", d);
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
        assert!(r.code == $code, "Retcodes are not equal.");
    };
}

#[macro_export]
macro_rules! test_check_result_sgn_and_code {
    ($func:ident, ($($arg:expr ),*), $val:expr, $tol:expr, $sgn:expr, $code:expr) => {
        let (r, sgn) = $func($($arg),*);
        let mut s:usize = 0;
        let mut f: f64 = 0.0;
        let mut d: f64 = 0.0;
        let ival :f64 = ($val as f64);

        if r.val.is_nan() ||ival.is_nan() {
            s = if r.val.is_nan() != ival.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else if r.val.is_infinite() || ival.is_infinite() {
            s = if r.val.is_infinite() != ival.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else {
            f = test_sf_frac_diff(ival, r.val);
            d = (ival - r.val).abs();

            if d > 2.5 * r.err {
                s |= TestResult::InCons as usize;
            }
            if r.err < 0.0 {
                s |= TestResult::ErrNeg as usize;
            }
            if r.err.is_infinite() {
                s |= TestResult::ErrBad as usize;
            }

            if f > $tol {
                s |= TestResult::TolBad as usize;
            }
        }

        if s != 0 {
            println!("expected: {:?}", ival);
            println!(
                "obtained: {:?} +/- {:?} (rel={:?})",
                r.val,
                r.err,
                r.err / (r.val + r.err).abs()
            );
            println!("fracdiff: {:?}", f);
            println!("tolerance: {:?}", $tol);
            println!("absdiff: {:?}", d);
            println!("call: {:?}({:?})",stringify!($func),($($arg),*));
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
        assert!(r.code == $code, "Retcodes are not equal.");
        assert!(sgn == $sgn, "Signs are not equal.");
    };

    (($res:expr, $sgn:expr), $val:expr, $tol:expr, $sgn_exp:expr, $code:expr) =>{
        let r = $res;
        let mut s:usize = 0;
        let mut f: f64 = 0.0;
        let mut d: f64 = 0.0;
        let ival :f64 = ($val as f64);

        if r.val.is_nan() ||ival.is_nan() {
            s = if r.val.is_nan() != ival.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else if r.val.is_infinite() || ival.is_infinite() {
            s = if r.val.is_infinite() != ival.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else {
            f = test_sf_frac_diff(ival, r.val);
            d = (ival - r.val).abs();

            if d > 2.0 * r.err {
                s |= TestResult::InCons as usize;
            }
            if r.err < 0.0 {
                s |= TestResult::ErrNeg as usize;
            }
            if r.err.is_infinite() {
                s |= TestResult::ErrBad as usize;
            }

            if f > $tol {
                s |= TestResult::TolBad as usize;
            }
        }

        if s != 0 {
            println!("expected: {:?}", ival);
            println!(
                "obtained: {:?} +/- {:?} (rel={:?})",
                r.val,
                r.err,
                r.err / (r.val + r.err).abs()
            );
            println!("fracdiff: {:?}", f);
            println!("tolerance: {:?}", $tol);
            println!("absdiff: {:?}", d);
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
        assert!(r.code == $code, "Retcodes are not equal.");
        assert!($sgn == $sgn_exp, "Signs are not equal.");
    };
}

#[macro_export]
macro_rules! test_check_result_and_code_e10 {
    ($func:ident, ($($arg:expr ),*), $val:expr, $e10:expr, $tol:expr, $code:expr) => {
        let r = $func($($arg),*);
        let mut s:usize = 0;
        let mut f: f64 = 0.0;
        let mut d: f64 = 0.0;
        let ival :f64 = ($val as f64);

        if r.val.is_nan() ||ival.is_nan() {
            s = if r.val.is_nan() != ival.is_nan() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else if r.val.is_infinite() || ival.is_infinite() {
            s = if r.val.is_infinite() != ival.is_infinite() {
            TestResult::InCons as usize
        } else {
            s
        };
        } else {
        f = test_sf_frac_diff(ival, r.val);
        d = (ival - r.val).abs();

        if d > 2.0 * r.err {
            s |= TestResult::InCons as usize;
        }
        if r.err < 0.0 {
            s |= TestResult::ErrNeg as usize;
        }
        if r.err.is_infinite() {
            s |= TestResult::ErrBad as usize;
        }

        if f > $tol {
            s |= TestResult::TolBad as usize;
        }
    }

    if s != 0 {
        println!("expected: {:?}", ival);
        println!(
            "obtained: {:?} +/- {:?} (rel={:?})",
            r.val,
            r.err,
            r.err / (r.val + r.err).abs()
        );
        println!("fracdiff: {:?}", f);
        println!("tolerance: {:?}", $tol);
        println!("absdiff: {:?}", d);
        println!("call: {:?}({:?})",stringify!($func),($($arg),*));
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
    assert!(r.code == $code, "Retcodes are not equal.");

    let s = if r.e10 != $e10 {
        TestResult::ExpBad as usize
    } else {
        0
    };

    if s != 0 {
        println!("Expected exponent: 10^{:?}", $e10);
        println!("Obtained exponent: 10^{:?}", r.e10);
    }

    assert!(
        s & TestResult::ExpBad as usize == 0,
        "Exponent is incorrect."
    );
    };
}
