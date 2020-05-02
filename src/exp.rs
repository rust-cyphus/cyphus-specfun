pub(crate) mod core;

#[cfg(test)]
mod test {
    use super::core::*;
    use crate::consts::{LN_DBL_MAX, SQRT_DLB_EPS};
    #[macro_use]
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_check_result_and_code_e10;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;
    #[test]
    fn test_exp_e() {
        test_check_result_and_code!(
            exp_e,
            (-10.0),
            (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(exp_e, (10.0), (10.0_f64).exp(), TOL0, SpecFunCode::Success);
    }
    #[test]
    fn test_exp_e10_e() {
        test_check_result_and_code_e10!(
            exp_e10_e,
            (1.0),
            std::f64::consts::E,
            0,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (2000.0),
            3.88118019428363725,
            868,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (100.0),
            2.688117141816135448412625551e43,
            0,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (1000.0),
            1.970071114017046993888879352,
            434,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (-100.0),
            3.720075976020835962959695803e-44,
            0,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (-1000.0),
            5.075958897549456765291809479,
            -435,
            TOL3,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_err_e() {
        test_check_result_and_code!(
            exp_err_e,
            (-10.0, TOL1),
            (-10.0_f64).exp(),
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_err_e,
            (10.0, TOL1),
            (10.0_f64).exp(),
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_err_e10_e() {
        test_check_result_and_code_e10!(
            exp_err_e10_e,
            (1.0, SQRT_TOL0),
            std::f64::consts::E,
            0,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_err_e10_e,
            (2000.0, 1e-10),
            3.88118019428363725,
            868,
            1e-7,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_mult_e() {
        let x = 0.8 * LN_DBL_MAX;
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, 1.0e-06),
            1.0e-06 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, 2.0),
            2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, -2.0),
            -2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (10.0, 1.0e-06),
            1.0e-06 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (10.0, -2.0),
            -2.0 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.00001),
            1.00001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.000001),
            1.000001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.000000001),
            1.000000001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 100.0),
            100.0 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.0e+20),
            1.0e+20 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, (-x).exp() * std::f64::consts::LN_2.exp()),
            2.0,
            TOL4,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exp_mult_err_e() {
        let x = 0.8 * LN_DBL_MAX;
        test_check_result_and_code!(
            exp_mult_err_e,
            (-10.0, SQRT_TOL0, 2.0, SQRT_TOL0),
            2.0 * (-10.0_f64).exp(),
            SQRT_TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_err_e,
            (
                x,
                SQRT_TOL0 * x,
                (-x).exp() * std::f64::consts::LN_2.exp(),
                SQRT_TOL0 * (-x).exp() * std::f64::consts::LN_2.exp()
            ),
            2.0,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exp_mult_e10_e() {
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (1.0, 1.0),
            std::f64::consts::E,
            0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (1000.0, 1e200),
            1.970071114017046993888879352,
            634,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (10000.0, 1.0),
            8.806818225662921587261496007,
            4342,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (100.0, 1.0),
            2.688117141816135448412625551e43,
            0,
            TOL2,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exp_mult_err_e10_e() {
        test_check_result_and_code_e10!(
            exp_mult_err_e10_e,
            (1.0, TOL0, 1.0, TOL0),
            std::f64::consts::E,
            0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_err_e10_e,
            (1000.0, 1e-12, 1e200, 1e190),
            1.9700711140165661,
            634,
            TOL3,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_expm1_e() {
        test_check_result_and_code!(
            expm1_e,
            (-10.0),
            (-10.0_f64).exp() - 1.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (-0.001),
            -0.00099950016662500845,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (-1.0e-8),
            -1.0e-08 + 0.5e-16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (1.0e-8),
            1.0e-08 + 0.5e-16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (0.001),
            0.0010005001667083417,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (10.0),
            (10.0_f64).exp() - 1.0,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exprel_e() {
        test_check_result_and_code!(
            exprel_e,
            (-10.0),
            0.0999954600070237515,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (-0.001),
            0.9995001666250084,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (-1.0e-8),
            1.0 - 0.5e-08,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (1.0e-8),
            1.0 + 0.5e-08,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (0.001),
            1.0005001667083417,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (10.0),
            2202.5465794806716517,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exprel_2_e() {
        test_check_result_and_code!(
            exprel_2_e,
            (-10.0),
            0.18000090799859524970,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (-0.001),
            0.9996667499833361107,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (-1.0e-8),
            0.9999999966666666750,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (1.0e-8),
            1.0000000033333333417,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (0.001),
            1.0003334166833361115,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (10.0),
            440.3093158961343303,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_exprel_n_e() {
        test_check_result_and_code!(
            exprel_n_e,
            (3, -1000.0),
            0.00299400600000000000,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -100.0),
            0.02940600000000000000,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -10.0),
            0.24599972760042142509,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -3.0),
            0.5444917625849191238,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -0.001),
            0.9997500499916678570,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -1.0e-8),
            0.9999999975000000050,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 1.0e-8),
            1.0000000025000000050,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 0.001),
            1.0002500500083345240,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 3.0),
            2.5745637607083706091,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 3.1),
            2.6772417068460206247,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 10.0),
            131.79279476884029910,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 100.0),
            1.6128702850896812690e+38,
            TOL2,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            exprel_n_e,
            (50, -1000.0),
            0.04766231609253975959,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -100.0),
            0.3348247572345889317,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -10.0),
            0.8356287051853286482,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -3.0),
            0.9443881609152163615,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -1.0),
            0.980762245565660617,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -1.0e-8),
            1.0 - 1.0e-8 / 51.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 1.0e-8),
            1.0 + 1.0e-8 / 51.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 1.0),
            1.01999216583666790,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 3.0),
            1.0624205757460368307,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 48.0),
            7.499573876877194416,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 50.1),
            9.311803306230992272,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 100.0),
            8.175664432485807634e+07,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 500.0),
            4.806352370663185330e+146,
            TOL3,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            exprel_n_e,
            (500, -1000.0),
            0.3334815803127619256,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -100.0),
            0.8335646217536183909,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -10.0),
            0.9804297803131823066,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -3.0),
            0.9940475488850672997,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -1.0),
            0.9980079602383488808,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -1.0e-8),
            1.0 - 1.0e-8 / 501.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1.0e-8),
            1.0 + 1.0e-8 / 501.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1.0),
            1.0019999920160634252,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 3.0),
            1.0060240236632444934,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 48.0),
            1.1059355517981272174,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 100.0),
            1.2492221464878287204,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 500.0),
            28.363019877927630858,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1000.0),
            2.4037563160335300322e+68,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1600.0),
            7.899293535320607403e+226,
            TOL4,
            SpecFunCode::Success
        );
    }
}
