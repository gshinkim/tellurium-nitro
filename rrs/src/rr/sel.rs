//! roadrunner selection strings: `time`, `S` (amount, or any symbol's value), `[S]`, `init(x)`,
//! `init([S])`, `x'`, `[S]'`, `eigen(S)`, `eigenReal(S)`, `eigenImag(S)`, `cc(a, b)`, `ucc(a, b)`,
//! `ec(r, x)`, `uec(r, x)`, `stoich(S, R)`. Parsing is syntax only; ids are resolved by the caller.

#[derive(Debug, Clone, PartialEq)]
pub enum Sel {
    /// Value of a symbol (species: amount; `time`: the model time unless a model id is named so).
    Id(String),
    /// `[S]`: species concentration.
    Conc(String),
    /// `init(x)` (false) / `init([S])` (true).
    Init(String, bool),
    /// `x'` (false) / `[S]'` (true): rate of change.
    Rate(String, bool),
    /// eigenvalue of the species' row: 0 = `eigen` (real part as a number), 1 = real, 2 = imaginary.
    Eigen(String, u8),
    /// `cc(a, b)` (true) / `ucc(a, b)`.
    Cc(String, String, bool),
    /// `ec(r, x)` (true) / `uec(r, x)`.
    Ec(String, String, bool),
    Stoich(String, String),
}

fn id(s: &str) -> Result<String, String> {
    let s = s.trim();
    if s.is_empty() || s.contains(['(', ')', '[', ']', ',', '\'', ' ']) { return Err(format!("invalid symbol '{s}'")); }
    Ok(s.to_string())
}

/// `[x]` -> (x, true), `x` -> (x, false).
fn br(s: &str) -> Result<(String, bool), String> {
    let s = s.trim();
    match s.strip_prefix('[').and_then(|x| x.strip_suffix(']')) {
        Some(x) => Ok((id(x)?, true)),
        None => Ok((id(s)?, false)),
    }
}

pub fn parse(s: &str) -> Result<Sel, String> {
    let t = s.trim();
    let bad = || format!("invalid selection '{s}'");
    if let Some(x) = t.strip_suffix('\'') {
        let (x, c) = br(x).map_err(|_| bad())?;
        return Ok(Sel::Rate(x, c));
    }
    if let Some(open) = t.find('(') {
        let (f, rest) = (t[..open].trim(), &t[open + 1..]);
        let args: Vec<&str> = rest.strip_suffix(')').ok_or_else(bad)?.split(',').collect();
        let two = |c: fn(String, String) -> Sel| -> Result<Sel, String> {
            if args.len() != 2 { return Err(bad()); }
            Ok(c(id(args[0]).map_err(|_| bad())?, id(args[1]).map_err(|_| bad())?))
        };
        let one = || if args.len() == 1 { id(args[0]).map_err(|_| bad()) } else { Err(bad()) };
        return match f {
            "init" if args.len() == 1 => { let (x, c) = br(args[0]).map_err(|_| bad())?; Ok(Sel::Init(x, c)) }
            "eigen" => Ok(Sel::Eigen(one()?, 0)),
            "eigenReal" => Ok(Sel::Eigen(one()?, 1)),
            "eigenImag" => Ok(Sel::Eigen(one()?, 2)),
            "cc" => two(|a, b| Sel::Cc(a, b, true)),
            "ucc" => two(|a, b| Sel::Cc(a, b, false)),
            "ec" => two(|a, b| Sel::Ec(a, b, true)),
            "uec" => two(|a, b| Sel::Ec(a, b, false)),
            "stoich" => two(Sel::Stoich),
            _ => Err(bad()),
        };
    }
    let (x, c) = br(t).map_err(|_| bad())?;
    Ok(if c { Sel::Conc(x) } else { Sel::Id(x) })
}

#[cfg(test)]
mod tests {
    use super::*;
    use Sel::*;

    #[test]
    fn selections() {
        let s = |x: &str| x.to_string();
        assert_eq!(parse("time"), Ok(Id(s("time"))));
        assert_eq!(parse(" S1 "), Ok(Id(s("S1"))));
        assert_eq!(parse("[S1]"), Ok(Conc(s("S1"))));
        assert_eq!(parse("init(S1)"), Ok(Init(s("S1"), false)));
        assert_eq!(parse("init([S1])"), Ok(Init(s("S1"), true)));
        assert_eq!(parse("init(k1)"), Ok(Init(s("k1"), false)));
        assert_eq!(parse("S1'"), Ok(Rate(s("S1"), false)));
        assert_eq!(parse("[S1]'"), Ok(Rate(s("S1"), true)));
        assert_eq!(parse("eigen(S1)"), Ok(Eigen(s("S1"), 0)));
        assert_eq!(parse("eigenReal(S1)"), Ok(Eigen(s("S1"), 1)));
        assert_eq!(parse("eigenImag(S1)"), Ok(Eigen(s("S1"), 2)));
        assert_eq!(parse("cc(J1, k1)"), Ok(Cc(s("J1"), s("k1"), true)));
        assert_eq!(parse("ucc(S1,k1)"), Ok(Cc(s("S1"), s("k1"), false)));
        assert_eq!(parse("ec(J1, S1)"), Ok(Ec(s("J1"), s("S1"), true)));
        assert_eq!(parse("uec(J1,k1)"), Ok(Ec(s("J1"), s("k1"), false)));
        assert_eq!(parse("stoich(S1, J1)"), Ok(Stoich(s("S1"), s("J1"))));
        for bad in ["", "[S1", "S1]", "init()", "init(S1", "cc(J1)", "foo(S1)", "[a b]", "cc(J1,k1,x)", "eigen(S1,S2)", "''"] {
            assert!(parse(bad).is_err(), "{bad} should not parse");
        }
    }
}
