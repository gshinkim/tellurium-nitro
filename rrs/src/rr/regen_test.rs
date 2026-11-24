//! regenerate() as libroadrunner's regenerateModel (values checked against libroadrunner 2.11).
use super::RoadRunner;

fn load(src: &str) -> RoadRunner {
    let sbml = antimony::sbml::write_sbml(&antimony::parse(src).document);
    let mut r = RoadRunner::new();
    r.load_sbml(&sbml, None).unwrap();
    r
}

#[test]
fn regenerate_carries_current_values() {
    let mut r = load("J1: S1 -> S2; k*S1; S1 = 10; S2 = 0; k = 0.5");
    (r.end, r.points) = (1.0, 3);
    r.simulate().unwrap();
    r.set_value("k", 0.9).unwrap();
    let s1 = r.get_value("S1").unwrap();
    // an edit that changes k's initial value: current values and time still carry over
    r.edit(|s| { assert!(s.contains("value=\"0.5\"")); *s = s.replacen("value=\"0.5\"", "value=\"0.2\"", 1); Ok(()) }, true).unwrap();
    assert_eq!((r.get_value("k").unwrap(), r.get_value("S1").unwrap(), r.get_value("time").unwrap()), (0.9, s1, 1.0));
    assert_eq!(r.get_value("init(k)").unwrap(), 0.2);
    r.reset().unwrap();
    assert_eq!((r.get_value("S1").unwrap(), r.get_value("k").unwrap()), (10.0, 0.9));
}
