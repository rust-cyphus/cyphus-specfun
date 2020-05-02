use crate::result::{SpecFunCode, SpecFunResult};

use crate::zeta::data::*;

pub fn hzeta_e(s: f64, q: f64) -> SpecFunResult<f64> {
    if s <= 1.0 || q <= 0.0 {
        let result = SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("hzeta_e", &[s, q]);
        result
    } else {
        let max_bits = 54.0;
        let ln_term0 = -s * q.ln();

        if ln_term0 < std::f64::MIN_POSITIVE.ln() + 1.0 {
            let result = SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::UnderflowErr,
            };
            result.issue_warning("hzeta_e", &[s, q]);
            result
        } else if ln_term0 > std::f64::MAX.ln() - 1.0 {
            let result = SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::OverflowErr,
            };
            result.issue_warning("hzeta_e", &[s, q]);
            result
        } else if (s > max_bits && q < 1.0) || (s > 0.5 * max_bits && q < 0.25) {
            let val = q.powf(-s);
            let err = 2.0 * std::f64::EPSILON * val.abs();
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        } else if s > 0.5 * max_bits && q < 1.0 {
            let p1 = q.powf(-s);
            let p2 = (q / (1.0 + q)).powf(s);
            let p3 = (q / (2.0 + q)).powf(s);
            let val = p1 * (1.0 + p2 + p3);
            let err = std::f64::EPSILON * (0.5 * s + 2.0) * val.abs();
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        } else {
            // Euler-Maclaurin summation formula
            // [Moshier, p. 400, with several typo corrections]
            let jmax = 12;
            let kmax = 10;
            let pmax = (kmax as f64 + q).powf(-s);
            let mut scp = s;
            let mut pcp = pmax / (kmax as f64 + q);
            let mut ans = pmax * ((kmax as f64 + q) / (s - 1.0) + 0.5);

            ans += (0..kmax).fold(0.0, |acc, k| acc + (k as f64 + q).powf(-s));

            for j in 0..(jmax + 1) {
                let delta = HZETA_C[j + 1] * scp * pcp;
                ans += delta;
                if (delta / ans).abs() < 0.5 * std::f64::EPSILON {
                    break;
                }
                scp *= (s + (2 * j) as f64 + 1.0) * (s + (2 * j) as f64 + 2.0);
                pcp /= (kmax as f64 + q) * (kmax as f64 + q);
            }
            let err = 2.0 * (jmax as f64 + 1.0) * std::f64::EPSILON * ans.abs();
            SpecFunResult {
                val: ans,
                err,
                code: SpecFunCode::Success,
            }
        }
    }
}
