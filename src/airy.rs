mod airy_core;
mod data;
use crate::result::SpecFunResult;

use airy_core::*;

pub trait Airy {
    /// Compute the Airy function Ai(x) along with an error estimate.
    fn airy_ai_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Ai(x) along with an error estimate,
    /// scaled by S(x) = exp(2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_ai_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Ai(x).
    fn airy_ai(&self) -> f64;
    /// Compute the Airy function Ai(x) scaled by S(x) = exp(2x^{3/2}/3)
    /// for x > 0 and 1 for x < 0.
    fn airy_ai_scaled(&self) -> f64;
    /// Compute the Airy function Bi(x) along with an error estimate.
    fn airy_bi_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Bi(x) along with an error estimate,
    /// scaled by S(x) = exp(-2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_bi_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Bi(x).
    fn airy_bi(&self) -> f64;
    /// Compute the Airy function Bi(x) scaled by S(x) = exp(-2x^{3/2}/3)
    /// for x > 0 and 1 for x < 0.
    fn airy_bi_scaled(&self) -> f64;
}

impl Airy for f64 {
    fn airy_ai_e(&self) -> SpecFunResult<f64> {
        airy_ai_e(*self)
    }
    fn airy_ai_scaled_e(&self) -> SpecFunResult<f64> {
        airy_ai_scaled_e(*self)
    }
    fn airy_ai(&self) -> f64 {
        airy_ai_e(*self).val
    }
    fn airy_ai_scaled(&self) -> f64 {
        airy_ai_scaled_e(*self).val
    }
    fn airy_bi_e(&self) -> SpecFunResult<f64> {
        airy_bi_e(*self)
    }
    fn airy_bi_scaled_e(&self) -> SpecFunResult<f64> {
        airy_bi_scaled_e(*self)
    }
    fn airy_bi(&self) -> f64 {
        airy_bi_e(*self).val
    }
    fn airy_bi_scaled(&self) -> f64 {
        airy_bi_scaled_e(*self).val
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;
    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * crate::consts::SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_airy_ai_e() {
        test_check_result_and_code!(
            airy_ai_e,
            (-500.0),
            0.0725901201040411396,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-5.0),
            0.3507610090241142,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-0.3000000000000094),
            0.4309030952855831,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (0.6999999999999907),
            0.1891624003981519,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (1.649999999999991),
            0.0583105861872088521,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (2.54999999999999),
            0.01446149513295428,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (3.499999999999987),
            0.002584098786989702,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (5.39999999999998),
            4.272986169411866e-05,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_scaled_e() {
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (-5.0),
            0.3507610090241142,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (0.6999999999999907),
            0.2795125667681217,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (1.649999999999991),
            0.2395493001442741,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (2.54999999999999),
            0.2183658595899388,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (3.499999999999987),
            0.2032920808163519,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (5.39999999999998),
            0.1836050093282229,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_e() {
        test_check_result_and_code!(
            airy_bi_e,
            (-500.0),
            -0.094688570132991028,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (-5.0),
            -0.1383691349016005,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (0.6999999999999907),
            0.9733286558781599,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (1.649999999999991),
            2.196407956850028,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (2.54999999999999),
            6.973628612493443,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (3.499999999999987),
            33.05550675461069,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (5.39999999999998),
            1604.476078241272,
            TOL1,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_scaled_e() {
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (-5.0),
            -0.1383691349016005,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (0.6999999999999907),
            0.6587080754582302,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (1.649999999999991),
            0.5346449995597539,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (2.54999999999999),
            0.461835455542297,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (3.499999999999987),
            0.4201771882353061,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (5.39999999999998),
            0.3734050675720473,
            TOL0,
            SpecFunCode::Success
        );
    }
}
