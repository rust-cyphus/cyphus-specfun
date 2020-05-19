pub(crate) mod besselj;
pub(crate) mod data;
mod helpers;

#[cfg(test)]
mod test {
    use crate::bessel::besselj::*;
    use crate::consts::SQRT_DLB_EPS;
    #[macro_use]
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_besselj0_e() {
        test_check_result_and_code!(
            besselj0_e,
            (0.1),
            0.99750156206604003230,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (2.0),
            0.22389077914123566805,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (100.0),
            0.019985850304223122424,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (1.0e+10),
            2.1755917502468917269e-06,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_besselj1_e() {
        test_check_result_and_code!(
            besselj1_e,
            (0.1),
            0.04993752603624199756,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (2.0),
            0.57672480775687338720,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (100.0),
            -0.07714535201411215803,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (1.0e+10),
            -7.676508175684157103e-06,
            TOL4,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_besseljn_e() {
        test_check_result_and_code!(
            besseljn_e,
            (4, 0.1),
            2.6028648545684032338e-07,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (5, 2.0),
            0.007039629755871685484,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (10, 20.0),
            0.18648255802394508321,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (100, 100.0),
            0.09636667329586155967,
            TOL0,
            SpecFunCode::Success
        );
    }
}
