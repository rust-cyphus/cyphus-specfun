use crate::cheb::ChebSeries;
use crate::result::{SpecFunCode, SpecFunResult};
use lazy_static::lazy_static;

lazy_static! {
    static ref ADEB1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.4006597190381410194,
            0.1937213042189360089,
            -0.62329124554895770e-02,
            0.3511174770206480e-03,
            -0.228222466701231e-04,
            0.15805467875030e-05,
            -0.1135378197072e-06,
            0.83583361188e-08,
            -0.6264424787e-09,
            0.476033489e-10,
            -0.36574154e-11,
            0.2835431e-12,
            -0.221473e-13,
            0.17409e-14,
            -0.1376e-15,
            0.109e-16,
            -0.9e-18,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB2_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.5943810232570770282,
            0.2863357204530719834,
            -0.102062656158046713e-01,
            0.6049109775346844e-03,
            -0.405257658950210e-04,
            0.28633826328811e-05,
            -0.2086394303065e-06,
            0.155237875826e-07,
            -0.11731280087e-08,
            0.897358589e-10,
            -0.69317614e-11,
            0.5398057e-12,
            -0.423241e-13,
            0.33378e-14,
            -0.2645e-15,
            0.211e-16,
            -0.17e-17,
            0.1e-18,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB3_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.707737068327440945,
            0.340068135211091751,
            -0.12945150184440869e-01,
            0.7963755380173816e-03,
            -0.546360009590824e-04,
            0.39243019598805e-05,
            -0.2894032823539e-06,
            0.217317613962e-07,
            -0.16542099950e-08,
            0.1272796189e-09,
            -0.987963460e-11,
            0.7725074e-12,
            -0.607797e-13,
            0.48076e-14,
            -0.3820e-15,
            0.305e-16,
            -0.24e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB4_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.781869415020523460,
            0.374976783526892863,
            -0.14940907399031583e-01,
            0.945679811437042e-03,
            -0.66132916138933e-04,
            0.4815632982144e-05,
            -0.3588083958759e-06,
            0.271601187416e-07,
            -0.20807099122e-08,
            0.1609383869e-09,
            -0.125470979e-10,
            0.9847265e-12,
            -0.777237e-13,
            0.61648e-14,
            -0.4911e-15,
            0.393e-16,
            -0.32e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB5_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.8340269546834530149,
            0.3994098857106266445,
            -0.164566764773099646e-1,
            0.10652138340664541e-2,
            -0.756730374875418e-4,
            0.55745985240273e-5,
            -0.4190692330918e-6,
            0.319456143678e-7,
            -0.24613318171e-8,
            0.1912801633e-9,
            -0.149720049e-10,
            0.11790312e-11,
            -0.933329e-13,
            0.74218e-14,
            -0.5925e-15,
            0.475e-16,
            -0.39e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB6_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.8726727134130122113,
            0.4174375352339027746,
            -0.176453849354067873e-1,
            0.11629852733494556e-2,
            -0.837118027357117e-4,
            0.62283611596189e-5,
            -0.4718644465636e-6,
            0.361950397806e-7,
            -0.28030368010e-8,
            0.2187681983e-9,
            -0.171857387e-10,
            0.13575809e-11,
            -0.1077580e-12,
            0.85893e-14,
            -0.6872e-15,
            0.552e-16,
            -0.44e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
}

pub fn debye1_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 1.64493406684822644;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 0.25 * x + x * x / 36.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB1_CHEB).eval(t);
        result.val = c.val - 0.25 * x;
        result.err = c.err + 0.25 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for i in (1..(nexp + 1)).rev() {
            sum *= ex;
            sum += 1f64.mul_add(xk.recip(), 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / x - sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        result.val = (val_infinity - (-x).exp() * (x + 1.0)) / x;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye2_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 4.80822761263837714;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - x / 3.0 + x * x / 24.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB2_CHEB).eval(t);
        result.val = c.val - x / 3.0;
        result.err = c.err + x / 3.0 * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for i in (1..(nexp + 1)).rev() {
            let xk_inv = xk.recip();
            sum *= ex;
            sum += 2f64.mul_add(xk_inv, 2.0).mul_add(xk_inv, 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x) - 2.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let sum = 2.0 + 2.0 * x + x2;
        result.val = (val_infinity - 2.0 * sum * (-x).exp()) / x2;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye3_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 19.4818182068004875;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 3.0 * x / 8.0 + x * x / 20.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB3_CHEB).eval(t);
        result.val = c.val - 0.375 * x;
        result.err = c.err + 0.375 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 6f64
                .mul_add(xk_inv, 6.0)
                .mul_add(xk_inv, 3.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x) - 3.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x3 = x * x * x;
        let sum = 6.0 + 6.0 * x + 3.0 * x * x + x3;
        result.val = (val_infinity - 3.0 * sum * (-x).exp()) / x3;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye4_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 99.5450644937635129;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * std::f64::consts::SQRT_2 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 2.0 * x / 5.0 + x * x / 18.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB4_CHEB).eval(t);
        result.val = c.val - 2.0 / 5.0 * x;
        result.err = c.err + 2.0 / 5.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 24f64
                .mul_add(xk_inv, 24.0)
                .mul_add(xk_inv, 12.0)
                .mul_add(xk_inv, 4.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x) - 4.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let sum = 24.0 + 24.0 * x + 12.0 * x2 + 4.0 * x2 * x + x4;
        result.val = (val_infinity - 4.0 * sum * (-x).exp()) / x4;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye5_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 610.405837190669483828710757875;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 5.0 * x / 12.0 + 5.0 * x * x / 84.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB5_CHEB).eval(t);
        result.val = c.val - 5.0 / 12.0 * x;
        result.err = c.err + 5.0 / 12.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 120f64
                .mul_add(xk_inv, 120.0)
                .mul_add(xk_inv, 60.0)
                .mul_add(xk_inv, 20.0)
                .mul_add(xk_inv, 5.0)
                .mul_add(xk_inv, 1.0)
                / rk;

            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x * x) - 5.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x5 = x4 * x;
        let sum = 120.0 + 120.0 * x + 60.0 * x2 + 20.0 * x2 * x + 5.0 * x4 + x5;
        result.val = (val_infinity - 5.0 * sum * (-x).exp()) / x5;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = ((((val_infinity / x) / x) / x) / x) / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye6_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 4356.06887828990661194792541535;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 3.0 * x / 7.0 + x * x / 16.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB6_CHEB).eval(t);
        result.val = c.val - 3.0 / 7.0 * x;
        result.err = c.err + 3.0 / 7.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 720f64
                .mul_add(xk_inv, 720.0)
                .mul_add(xk_inv, 360.0)
                .mul_add(xk_inv, 120.0)
                .mul_add(xk_inv, 30.0)
                .mul_add(xk_inv, 6.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x * x * x) - 6.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x4 * x2;
        let sum = 720.0 + 720.0 * x + 360.0 * x2 + 120.0 * x2 * x + 30.0 * x4 + 6.0 * x4 * x + x6;
        result.val = (val_infinity - 6.0 * sum * (-x).exp()) / x6;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = (((((val_infinity / x) / x) / x) / x) / x) / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;

    #[test]
    fn test_debye1_e() {
        test_check_result_and_code!(
            debye1_e,
            (0.1),
            0.975277750004723276,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye1_e,
            (1.0),
            0.777504634112248239,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye1_e,
            (10.0),
            0.164443465679946027,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_debye2_e() {
        test_check_result_and_code!(
            debye2_e,
            (0.1),
            0.967083287045302664,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye2_e,
            (1.0),
            0.70787847562782924,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye2_e,
            (10.0),
            0.0479714980201218708,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_debye3_e() {
        test_check_result_and_code!(
            debye3_e,
            (0.1),
            0.962999940487211048,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye3_e,
            (1.0),
            0.674415564077814667,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye3_e,
            (10.0),
            0.0192957656903454886,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_debye4_e() {
        test_check_result_and_code!(
            debye4_e,
            (0.1),
            0.960555486124335944,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye4_e,
            (1.0),
            0.654874068886737049,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye4_e,
            (10.0),
            0.00967367556027115896,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_debye5_e() {
        test_check_result_and_code!(
            debye5_e,
            (0.1),
            0.95892849428310568745,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye5_e,
            (1.0),
            0.6421002580217790246,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye5_e,
            (10.0),
            0.005701535852992908538,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_debye6_e() {
        test_check_result_and_code!(
            debye6_e,
            (0.1),
            0.95776777382605465878,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye6_e,
            (1.0),
            0.63311142583495107588,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye6_e,
            (10.0),
            3.7938493294615955279e-3,
            TOL0,
            SpecFunCode::Success
        );
    }
}
