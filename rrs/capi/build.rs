// Emits the C-style __DATE__ / __TIME__ strings (UTC) and the rustc version for getBuildDate & co.
use std::time::{SystemTime, UNIX_EPOCH};

fn main() {
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs() as i64;
    let (days, sod) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    // civil_from_days (H. Hinnant)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    const MON: [&str; 12] = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    println!("cargo:rustc-env=RRC_BUILD_DATE={} {:>2} {}", MON[(m - 1) as usize], d, y);
    println!("cargo:rustc-env=RRC_BUILD_TIME={:02}:{:02}:{:02}", sod / 3600, sod % 3600 / 60, sod % 60);
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    let ver = std::process::Command::new(rustc).arg("--version").output().ok()
        .and_then(|o| String::from_utf8(o.stdout).ok()).unwrap_or_else(|| "rustc".into());
    println!("cargo:rustc-env=RRC_RUSTC={}", ver.trim());
}
