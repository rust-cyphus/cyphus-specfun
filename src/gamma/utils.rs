use crate::{
    cheb::cheb_eval_e,
    consts::{LNPI, ROOT4_DBL_EPS},
    exp::core::{exp_err_e, exp_mult_err_e},
    result::{SpecFunCode, SpecFunResult},
    trig::angle_restrict_symm_e,
    zeta::Zeta,
};

use num::{complex::Complex, Float};

use super::data::*;

// ------------------
// Internal Functions
// ------------------

/// Compute the factorial of a number, returning an f64
pub(crate) fn factorial(n: usize) -> f64 {
    if n < FACT_TABLE.len() {
        FACT_TABLE[n]
    } else {
        std::f64::NAN
    }
}

/// Compute the double factorial of a number, returning an f64
pub(crate) fn double_factorial(n: usize) -> f64 {
    if n < DOUBLE_FACT_TABLE.len() {
        DOUBLE_FACT_TABLE[n]
    } else {
        std::f64::NAN
    }
}

/// Compute Log(Gamma(z)) using Lanczos method for complex numbers
pub(crate) fn lngamma_lanczos_complex_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut y = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(f64::EPSILON, f64::EPSILON),
        code: SpecFunCode::Success,
    };

    let zz = Complex::new(z.re - 1.0, z.im);
    let mut ag = Complex::new(LANCZOS_7_C[0], 0.0);
    for k in 1..9 {
        let r = zz.re + (k as f64);
        let i = zz.im;
        let a = LANCZOS_7_C[k] / (r * r + i * i);
        ag.re += a * r;
        ag.im -= a * i;
    }

    let log1 = (zz + 7.5).ln();
    let logag = ag.ln();

    y.val.re = (zz.re + 0.5) * log1.re - zz.im * log1.im - (zz.re + 7.5)
        + 0.918_938_533_204_672_741_8_f64
        + logag.re;
    y.val.im = zz.im * log1.re + (zz.re + 0.5) * log1.im - zz.im + logag.im;

    y.err.re = 4.0 * f64::EPSILON * y.val.re.abs();
    y.err.im = 4.0 * f64::EPSILON * y.val.im.abs();

    let res = angle_restrict_symm_e(y.val.im);
    y.val.im = res.val;
    y.err.im += res.err;
    y
}

/// Compute Log(Gamma(x)) using Lanczos method
pub(crate) fn lngamma_lanczos_e(x: f64) -> SpecFunResult<f64> {
    let xx = x - 1.0;
    let ag = LANCZOS_7_C
        .iter()
        .skip(1)
        .enumerate()
        .fold(LANCZOS_7_C[0], |acc, (k, val)| {
            acc + val / (xx + ((k + 1) as f64))
        });
    let term1 = (xx + 0.5) * ((xx + 7.5) * (-1.0).exp()).ln();
    let term2 = 0.918_938_533_204_672_741_8_f64 + ag.ln();

    let val = term1 + (term2 - 7.0);
    let err = f64::EPSILON * (2.0 * (term1.abs() + term2.abs() + 7.0) + val.abs());

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps / 2, as well as its sign
pub(crate) fn lngamma_sgn_0_e(eps: f64) -> (SpecFunResult<f64>, f64) {
    let c1 = - 0.077_215_664_901_532_860_61_f64;
    let c2 = - 0.010_944_004_672_027_444_61_f64;
    let c3 =   0.092_520_923_919_113_710_98_f64;
    let c4 = - 0.018_271_913_165_599_812_66_f64;
    let c5 =   0.018_004_931_096_854_797_90_f64;
    let c6 = - 0.006_850_885_378_723_806_85_f64;
    let c7 =   0.003_998_239_557_568_466_03_f64;
    let c8 = - 0.001_894_306_216_871_078_02_f64;
    let c9 =   0.000_974_732_378_045_132_21_f64;
    let c10 = -0.000_484_343_927_222_558_93_f64;
    let g6 = c6 + eps * (c7 + eps * (c8 + eps * (c9 + eps * c10)));
    let g = eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * (c5 + eps * g6)))));

    let gee = g + 1.0 / (1.0 + eps) + 0.5 * eps;
    let val = gee * eps.abs().recip();
    let err = 4.0 * f64::EPSILON * val.abs();
    (
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        },
        val.signum(),
    )
}

pub(crate) fn lngamma_sgn_sing_e(n: i32, eps: f64) -> (SpecFunResult<f64>, f64) {
    if eps.abs() < f64::EPSILON {
        (
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            },
            0.0,
        )
    } else if n == 1 {
        let c0 =  0.077_215_664_901_532_860_61_f64;
        let c1 =  0.088_159_669_573_560_305_21_f64;
        let c2 = -0.004_361_254_345_553_405_77_f64;
        let c3 =  0.013_910_658_820_046_406_89_f64;
        let c4 = -0.004_094_272_276_808_391_00_f64;
        let c5 =  0.002_756_613_101_915_415_84_f64;
        let c6 = -0.001_241_626_455_653_050_19_f64;
        let c7 =  0.000_652_679_761_218_027_83_f64;
        let c8 = -0.000_322_052_616_827_104_37_f64;
        let c9 =  0.000_162_291_310_395_454_56_f64;
        let g5 = c5 + eps * (c6 + eps * (c7 + eps * (c8 + eps * c9)));
        let g = eps * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * g5)))));

        // calculate eps gamma(-1+eps), a negative quantity
        let gam_e = g - 1.0 - 0.5 * eps * (1.0 + 3.0 * eps) / (1.0 - eps * eps);
        let val = (gam_e.abs() * eps.abs().recip()).ln();
        let err = 2.0 * f64::EPSILON * val.abs();
        (
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            },
            if eps > 0.0 { -1.0 } else { 1.0 },
        )
    } else {
        let cs1 = -1.644_934_066_848_226_436_5_f64;
        let cs2 =  0.811_742_425_283_353_643_6_f64;
        let cs3 = -0.190_751_824_122_084_213_7_f64;
        let cs4 =  0.026_147_847_817_654_800_5_f64;
        let cs5 = -0.002_346_081_035_455_823_6_f64;
        let e2 = eps * eps;
        let sin_ser = 1.0 + e2 * (cs1 + e2 * (cs2 + e2 * (cs3 + e2 * (cs4 + e2 * cs5))));

        // Calculate series for ln(gamma(1+n-eps))
        let aeps = eps.abs();

        let mut psi0 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi1 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi2 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi3 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi4 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi5 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };
        let mut psi6 = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };

        let c0 = lnfact_int_e(n as usize);
        psi0 = digamma_int_e((n + 1) as usize);
        psi1 = trigamma_int_e((n + 1) as usize);

        (
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            },
            0.0,
        )
    }
}

pub(crate) fn lnfact_int_e(n: usize) -> SpecFunResult<f64> {
    if n < FACT_TABLE.len() {
        let val = FACT_TABLE[n].ln();
        SpecFunResult {
            val: val,
            err: 2.0 * f64::EPSILON * val.abs(),
            code: SpecFunCode::Success,
        }
    } else {
        lngamma_e((n as f64) + 1.0)
    }
}

/// Compute log(Gamma(1+eps))/eps using the (2,2) pade
/// approximate plus a correction series
pub(crate) fn lngamma_1_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1 = -1.001_741_928_234_950_869_987_113_844_0_f64;
    let n2 =  1.736_483_920_992_287_982_328_054_173_3_f64;
    let d1 =  1.243_300_601_885_875_155_605_543_601_1_f64;
    let d2 =  5.045_627_410_027_401_015_248_959_751_4_f64;
    let num = (eps + n1) * (eps + n2);
    let den = (eps + d1) * (eps + d2);
    let pade = 2.081_626_518_866_269_247_488_021_031_8_f64 * num * den.recip();

    let c0 =  0.004_785_324_257_581_753_f64;
    let c1 = -0.011_924_570_836_454_41_f64;
    let c2 =  0.019_319_614_139_604_98_f64;
    let c3 = -0.025_940_273_987_250_20_f64;
    let c4 =  0.031_419_287_550_214_55_f64;
    let eps5 = eps * eps * eps * eps * eps;
    let corr = eps5 * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + c4 * eps))));

    // Error in case wanted for future use
    // let err = 2.0 * T::epsilon() * eps * (pade + corr);
    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();
    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute log(Gamma(2+eps))/eps using the (2,2) pade
/// approximate plus a correction series
pub(crate) fn lngamma_2_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1 = 1.000_895_834_786_669_227_164_446_568_f64;
    let n2 = 4.209_376_735_287_755_081_642_901_277_f64;
    let d1 = 2.618_851_904_903_217_274_682_578_255_f64;
    let d2 = 10.857_665_599_009_835_153_229_229_36_f64;
    let num = (eps + n1) * (eps + n2);
    let den = (eps + d1) * (eps + d2);
    let pade = 2.853_379_987_657_819_184_635_688_69_f64 * num / den;
    let c0 =   0.000_113_940_635_703_674_4_f64;
    let c1 = -0.000_136_543_526_979_253_3_f64;
    let c2 =   0.000_106_728_716_918_366_5_f64;
    let c3 = -0.000_069_327_180_093_128_2_f64;
    let c4 =   0.000_040_722_092_786_795_0_f64;
    let eps5 = eps * eps * eps * eps * eps;
    let corr = eps5 * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + c4 * eps))));

    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute the digamma function
fn digamma_int_e(n: usize) -> SpecFunResult<f64> {
    if n <= 0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("digamma_int_e", &[n as f64]);
        result
    } else if (n as usize) < PSI_TABLE.len(){
        let val = PSI_TABLE[n as usize];
        let err = std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        // Abramowitz+Stegun 6.3.18
        let nf = n as f64;
        let c2 = -1.0 / 12.0;
        let c3 = 1.0 / 120.0;
        let c4 = -1.0 / 252.0;
        let c5 = 1.0 / 240.0;
        let ni2 = (1.0 / nf) * (1.0 / nf);
        let ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2 * c5)));
        let val = nf.ln() - 0.5 / nf + ser;
        let mut err = std::f64::EPSILON * (nf.ln().abs() + (0.5 / nf).abs() + ser.abs());
        err += std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the trigamma function
fn trigamma_int_e(n: usize) -> SpecFunResult<f64> {
    if n <= 0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("trigamma_int_e", &[n as f64]);
        result
    } else if n as usize < PSI_1_TABLE.len() {
        let val = PSI_1_TABLE[n as usize];
        let err = std::f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        // Abramowitz+Stegun 6.4.12
        // double-precision for n > 100
        let nf = n as f64;
        let c0 = -1.0 / 30.0;
        let c1 = 1.0 / 42.0;
        let c2 = -1.0 / 30.0;
        let ni2 = (1.0 / nf) * (1.0 / nf);
        let ser = ni2 * ni2 * (c0 + ni2 * (c1 + c2 * ni2));
        let val = (1.0 + 0.5 / nf + 1.0 / (6.0 * nf * nf) + ser) / nf;
        let err = std::f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

pub(crate) fn gammastar_ser_e(x: f64) -> SpecFunResult<f64> {
    let y = 1.0 / (x * x);
    let c0 = 1.0 / 12.0;
    let c1 = -1.0 / 360.0;
    let c2 = 1.0 / 1_260.0;
    let c3 = -1.0 / 1_680.0;
    let c4 = 1.0 / 1_188.0;
    let c5 = -691.0 / 360_360.0;
    let c6 = 1.0 / 156.0;
    let c7 = -3617.0 / 122_400.0;
    let ser = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));

    let val = (ser / x).exp();
    let err = 2.0 * f64::EPSILON * val * (ser / x).max(1.0);

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute Gamma(x) for x >= 1/2
pub(crate) fn gamma_x_gt_half_e(x: f64) -> SpecFunResult<f64> {
    if (x - 0.5).abs() < f64::EPSILON {
        // Error term
        let val = 1.772_453_850_905_516_027_298_17_f64;
        let err = f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x <= 1.0 + (FACT_TABLE.len() as f64) - 1.0 && (x - x.floor()).abs() < f64::EPSILON {
        let n = x.floor() as usize;
        let val = FACT_TABLE[n - 1];
        let err = val * f64::EPSILON;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if (x - 1.0).abs() < 0.01 {
        let eps = x - 1.0;
        let c1 = 0.422_784_335_098_467_139_4_f64;
        let c2 = -0.010_944_004_672_027_444_61_f64;
        let c3 = 0.092_520_923_919_113_710_98_f64;
        let c4 = -0.018_271_913_165_599_812_664_f64;
        let c5 = 0.018_004_931_096_854_797_895_f64;
        let c6 = -0.006_850_885_378_723_806_846_f64;
        let c7 = 0.003_998_239_557_568_466_030_f64;
        let err = f64::EPSILON;
        let val = c7
            .mul_add(eps, c6)
            .mul_add(eps, c5)
            .mul_add(eps, c4)
            .mul_add(eps, c3)
            .mul_add(eps, c2)
            .mul_add(eps, c1)
            .mul_add(eps, x.recip());

        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if (x - 2.0).abs() < 0.01 {
        let eps = x - 2.0;
        let c1 = 0.422_784_335_098_467_139_4_f64;
        let c2 = 0.411_840_330_426_439_694_8_f64;
        let c3 = 0.081_576_919_247_086_266_38_f64;
        let c4 = 0.074_249_010_753_513_898_32_f64;
        let c5 = -0.000_266_982_068_745_014_768_32_f64;
        let c6 = 0.011_154_045_718_130_991_049_f64;
        let c7 = -0.002_852_645_821_155_340_816_f64;
        let c8 = 0.002_103_933_340_697_388_008_5_f64;

        let err = f64::EPSILON;
        let val = c8
            .mul_add(eps, c7)
            .mul_add(eps, c6)
            .mul_add(eps, c5)
            .mul_add(eps, c4)
            .mul_add(eps, c3)
            .mul_add(eps, c2)
            .mul_add(eps, c1)
            .mul_add(eps, 1.0);

        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 5.0 {
        // Exponentiating the logarithm is fine, as
        // long as the exponential is not so large
        // that it greatly amplifies the error.
        let mut lg = lngamma_lanczos_e(x);
        lg.val = lg.val.exp();
        lg.err = lg.val * (lg.err + 2.0 * f64::EPSILON);
        lg
    } else if x < 10.0 {
        // This is a sticky area. The logarithm
        // is too large and the gammastar series
        // is not good.
        let gamma_8 = 5040.0;
        let t = (2.0 * x - 15.0) / 5.0;
        let c = cheb_eval_e(t, &GAMMA_5_10_DATA, -1.0, 1.0);
        let val = c.val.exp() * gamma_8;
        let err = val * c.err + 2.0 * f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < GSL_SF_GAMMA_XMAX {
        // We do not want to exponentiate the logarithm
        // if x is large because of the inevitable
        // inflation of the error. So we carefully
        // use pow() and exp() with exact quantities.
        let p = x.powf(x * 0.5);
        let e = (-x).exp();
        let q = (p * e) * p;
        let pre = std::f64::consts::SQRT_2 * std::f64::consts::PI.sqrt() * q / x.sqrt();
        let mut gstar = gammastar_ser_e(x);

        gstar.val *= pre;
        gstar.err = (x + 2.5) * f64::EPSILON * gstar.val;

        gstar
    } else {
        let result = SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
            code: SpecFunCode::OverflowErr,
        };
        result.issue_warning("psi_x_e", &[x]);
        result
    }
}

/// Compute the digamma function for either positive or negative `x`.
pub(crate) fn psi_x_e(x: f64) -> SpecFunResult<f64> {
    let y = x.abs();

    if y < f64::EPSILON || (x + 1.0).abs() < f64::EPSILON || (x + 2.0).abs() < f64::EPSILON {
        let result = SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("psi_x_e", &[x]);
        result
    } else if y >= 2.0 {
        let t = 8.0 / (y * y) - 1.0;
        let result_c = cheb_eval_e(t, &APSICS_DATA, -1.0, 1.0);
        if x < 0.0 {
            let s = (x * std::f64::consts::PI).sin();
            let c = (x * std::f64::consts::PI).cos();
            if s.abs() < 2.0 * f64::MIN_POSITIVE.sqrt() {
                let result = SpecFunResult {
                    val: std::f64::NAN,
                    err: std::f64::NAN,
                    code: SpecFunCode::DomainErr,
                };
                result.issue_warning("psi_x_e", &[x]);
                result
            } else {
                let val = y.ln() - 0.5 / x + result_c.val - std::f64::consts::PI * c / s;
                let mut err = std::f64::consts::PI * x.abs() * f64::EPSILON / (s * s);
                err += result_c.err;
                err += f64::EPSILON * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            }
        } else {
            let val = y.ln() - 0.5 / x + result_c.val;
            let mut err = result_c.err;
            err += f64::EPSILON * val;
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        }
    } else if x < -1.0 {
        let v = x + 2.0;
        let t1 = 1.0 / x;
        let t2 = 1.0 / (x + 1.0);
        let t3 = 1.0 / v;
        let result_c = cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -(t1 + t2 + t3) + result_c.val;
        let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs() + (x / (t3 * t3)).abs());
        err += result_c.err;
        err += f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 0.0 {
        let v = x + 1.0;
        let t1 = 1.0 / x;
        let t2 = 1.0 / v;
        let result_c = cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -(t1 + t2) + result_c.val;
        let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs());
        err += result_c.err;
        err += f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 1.0 {
        let t1 = 1.0 / x;
        let result_c = cheb_eval_e(2.0 * x - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -t1 + result_c.val;
        let mut err = f64::EPSILON * t1;
        err += result_c.err;
        err += f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let v = x - 1.0;
        cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0)
    }
}

/// Compute psi(z) for large |z| in the right-half plane
/// ref: [Abramowitz + Stegun, 6.3.18]
pub(crate) fn psi_complex_asymp(z: Complex<f64>) -> Complex<f64> {
    // coefficients in the asymptotic expansion for large z;
    // let w = z^(-2) and write the expression in the form
    //
    //   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
    let c1 = -0.1;
    let c2 = 1.0 / 21.0;
    let c3 = -0.05;

    let zinv = z.inv();
    let w = zinv.powi(2);

    // Horner method evaluation of term in parentheses
    let mut sum: Complex<f64>;
    sum = w * (c3 / c2);
    sum = sum + 1.0;
    sum = sum * (c2 / c1);
    sum = sum * w;
    sum = sum + 1.0;
    sum = sum * c1;
    sum = sum * w;
    sum = sum + 1.0;

    // Correction added to log(z)
    let mut cs = sum * w;
    cs = cs * (-1.0 / 12.0);
    cs = cs + zinv * (-0.5);

    cs + z.ln()
}

/// Compute Psi(z) for complex z in the right-half plane
pub(crate) fn psi_complex_rhp(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut n_recurse: usize = 0;
    let mut res = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    if z.re.abs() < std::f64::EPSILON && z.im.abs() < std::f64::EPSILON {
        res.val = Complex::new(std::f64::NAN, std::f64::NAN);
        res.err = Complex::new(std::f64::NAN, std::f64::NAN);
        res.code = SpecFunCode::DomainErr;
        res.issue_warning("psi_complex_rhp", &[z]);
        return res;
    }

    // Compute the number of recurrences to apply
    if z.re < 20.0 && z.im.abs() < 20.0 {
        let sp = (20.0 + z.im).sqrt();
        let sn = (20.0 - z.im).sqrt();
        let rhs = sp * sn - z.re;
        if rhs > 0.0 {
            n_recurse = rhs.ceil() as usize;
        }
    }

    // Compute asymptotic at the large value z + n_recurse
    let mut a = psi_complex_asymp(z + n_recurse as f64);
    res.err = 2.0 * std::f64::EPSILON * Complex::new(a.re.abs(), a.im.abs());

    // Descend recursively, if necessary
    for i in (1..(n_recurse + 1)).rev() {
        let zn = z + (i as f64 - 1.0);
        let zn_inv = zn.inv();
        a = z - zn_inv;

        // Accumulate the error, to catch cancellations
        res.err.re += 2.0 * std::f64::EPSILON * zn_inv.re.abs();
        res.err.im += 2.0 * std::f64::EPSILON * zn_inv.im.abs();
    }

    res.val = a;

    res.err.re += 2.0 * std::f64::EPSILON * res.val.re.abs();
    res.err.im += 2.0 * std::f64::EPSILON * res.val.im.abs();

    res
}

/// Compute Psi^{(n)}(x) for x > 0
pub(crate) fn psi_n_xg0(n: usize, x: f64) -> SpecFunResult<f64> {
    if n == 0 {
        psi_x_e(x)
    } else {
        // Abramowitz + Stegun 6.4.10
        let hz = (n as f64 + 1.0).hzeta_e(x);
        let lnnf = lnfact_int_e(n);
        let mut result = exp_mult_err_e(lnnf.val, lnnf.err, hz.val, hz.err);

        if n % 2 == 0 {
            result.val *= -1.0;
        }

        result
    }
}

// ------------------
// External functions
// ------------------

/// ln(gamma(x)) where x is not a negative integer
///
// Uses real Lanczos method.
// Returns the real part of ln(gamma(x)) when x < 0,
// i.e. ln(|gamma(x)|).
pub(crate) fn lngamma_e(x: f64) -> SpecFunResult<f64> {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        result
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        result
    } else if x >= 0.5 {
        lngamma_lanczos_e(x)
    } else if x.abs() < f64::EPSILON {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_e", &[x]);
        result
    } else if x.abs() < 0.02 {
        lngamma_sgn_0_e(x).0
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try tp extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * x).sin();
        let abss = s.abs();
        if abss < f64::EPSILON {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("lngamma_e", &[x]);
            result
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::RoundoffErr,
                };
                result.issue_warning("lngamma_e", &[x]);
                result
            } else {
                let n = -(x - 0.5) as i32;
                let eps = x + n as f64;
                lngamma_sgn_sing_e(n, eps).0
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            result.val = LNPI - (abss.ln() + result.val);
            result.err = 2.0 * f64::EPSILON * result.val.abs() + result.err;
            result
        }
    } else {
        // |x| was too large to extract any fraction part
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::RoundoffErr,
        };
        result.issue_warning("lngamma_e", &[x]);
        result
    }
}

/// Compute ln(gamma(x)) where x is not a negative integer.
///
// Uses real Lanczos method. Determines
// the sign of gamma(x) as well as ln(|gamma(x)|) for x < 0.
// So gamma(x) = sgn * exp(lngamma_sgn_e(x).res).
pub(crate) fn lngamma_sgn_e(x: f64) -> (SpecFunResult<f64>, f64) {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        (result, 1.0)
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        (result, 1.0)
    } else if x >= 0.5 {
        (lngamma_lanczos_e(x), 1.0)
    } else if x.abs() < f64::EPSILON {
        let mut result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_sgn_e", &[x]);
        (result, 0.0)
    } else if x.abs() < 0.02 {
        lngamma_sgn_0_e(x)
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try to extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * z).sin();
        let abss = s.abs();
        if abss < f64::EPSILON {
            let mut result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("lngamma_sgn_e", &[x]);
            (result, 0.0)
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::RoundoffErr,
                };
                result.issue_warning("lngamma_sgn_e", &[x]);
                (result, 0.0)
            } else {
                let n = -(x - 0.5) as i32;
                let eps = x + n as f64;
                lngamma_sgn_sing_e(n, eps)
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            let sgn = if s > 0.0 { 1.0 } else { -1.0 };
            result.val = LNPI - (abss.ln() + result.val);
            result.err = 2.0 * f64::EPSILON * result.val.abs() + result.err;
            (result, sgn)
        }
    } else {
        // |x| was too large to extract any fraction part
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_sgn_e", &[x]);
        (result, 0.0)
    }
}
/// Gamma(x), where x is not a negative integer
///
/// Uses real Lanczos method.
pub(crate) fn gamma_e(x: f64) -> SpecFunResult<f64> {
    if x < 0.5 {
        let rint_x = (x + 0.5).floor() as i32;
        let f_x = x - rint_x as f64;
        let sgn_gamma = if rint_x % 2 == 0 { 1.0 } else { -1.0 };
        let sin_term = sgn_gamma * (std::f64::consts::PI * f_x).sin() / std::f64::consts::PI;

        if sin_term.abs() < f64::EPSILON {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("gamma_e", &[x]);
            result
        } else if x > -169.0 {
            let g = gamma_x_gt_half_e(1.0 - x);
            if sin_term.abs() * g.val * f64::MIN_POSITIVE < 1.0 {
                let val = 1.0 / (sin_term * g.val);
                let mut err = (g.err / g.val).abs() * val.abs();
                err += 2.0 * f64::EPSILON * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            } else {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::UnderflowErr,
                };
                result.issue_warning("gamma_e", &[x]);
                result
            }
        } else {
            // It is hard to control it here.
            // We can only exponentiate the
            // logarithm and eat the loss of
            // precision.
            let (lng, sgn) = lngamma_sgn_e(x);
            exp_mult_err_e(lng.val, lng.err, sgn, 0.0)
        }
    } else {
        gamma_x_gt_half_e(x)
    }
}

/// Regulated Gamma Function, x > 0
///
/// Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
///            = (1 + 1/(12x) + ...), x->Inf
pub(crate) fn gammastar_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("gammastar_e", &[x]);
        result
    } else if x < 0.5 {
        let lg = lngamma_e(x);
        let lx = x.ln();
        let c = 0.5 * (std::f64::consts::LN_2 + LNPI);
        let lnr_val = lg.val - (x - 0.5) * lx + x - c;
        let lnr_err = lg.err + 2.0 * f64::EPSILON * ((x + 0.5) * lx.abs() + c);
        exp_err_e(lnr_val, lnr_err)
    } else if x < 2.0 {
        let t = 4.0 / 3.0 * (x - 0.5) - 1.0;
        cheb_eval_e(t, &GSTAR_A_DATA, -1.0, 1.0)
    } else if x < 10.0 {
        let t = 0.25 * (x - 2.0) - 1.0;
        let c = cheb_eval_e(t, &gstar_b_data, -1.0, 1.0);
        let val = c.val / (x * x) + 1.0 + 1.0 / (12.0 * x);
        let mut err = c.err / (x * x);
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 1.0 / ROOT4_DBL_EPS {
        gammastar_ser_e(x)
    } else if x < 1.0 / f64::EPSILON {
        let xi = 1.0 / x;
        let val = 1.0
            + xi / 12.0 * (1.0 + xi / 24.0 * (1.0 - xi * (139.0 / 180.0 + 571.0 / 8640.0 * xi)));
        let err = 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        SpecFunResult {
            val: 1.0,
            err: x.recip(),
            code: SpecFunCode::Success,
        }
    }
}

/// Compute 1 / gamma(x)
///
/// Uses real Lanczos method.
pub(crate) fn gammainv_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 && (x - x.floor()).abs() < f64::EPSILON {
        SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        }
    } else if x < 0.5 {
        let (lng, sgn) = lngamma_sgn_e(x);
        if lng.val == f64::NAN {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        } else {
            exp_mult_err_e(-lng.val, lng.err, sgn, 0.0)
        }
    } else {
        let g = gamma_x_gt_half_e(x);
        let val = 1.0 / g.val;
        let mut err = (g.err / g.val).abs() * val.abs();
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fn() {
        assert!(true);
    }
}
