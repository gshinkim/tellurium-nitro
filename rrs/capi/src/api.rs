//! rrc_api.h on roadrunner_nitro::rr::RoadRunner. An RRHandle is a Box<Inst> (the engine plus the
//! C-side per-instance state). Semantics, ownership and failure values (false / NULL / -1) follow
//! libroadrunner's wrappers/C/rrc_api.cpp; errors go to the thread's getLastError, and no panic
//! crosses the FFI boundary.
use crate::types::*;
use crate::{has_error, last_error_ptr, set_error};
use libc::{c_char, c_int, c_uint};
use roadrunner_nitro::rr::{Matrix, RoadRunner, Solver, Value, INTEGRATORS, STEADY_STATE_SOLVERS};
use roadrunner_nitro::sbml::edit;
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::ptr::null_mut;
use std::sync::Mutex;

// ---------------------------------------------------------------- version / build info

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn version_ex() -> String {
    format!("{VERSION}, compiled with {} on {}, {}", env!("RRC_RUSTC"), env!("RRC_BUILD_DATE"), env!("RRC_BUILD_TIME"))
}

#[no_mangle]
pub extern "C" fn getAPIVersion() -> *mut c_char {
    text("1.0")
}

#[no_mangle]
pub extern "C" fn getCPPAPIVersion(_handle: RRHandle) -> *mut c_char {
    text(VERSION)
}

/// major*10000 + minor*100 + patch, as upstream.
#[no_mangle]
pub extern "C" fn getVersion() -> c_int {
    let p = |s: &str| s.parse::<c_int>().unwrap_or(0);
    p(env!("CARGO_PKG_VERSION_MAJOR")) * 10000 + p(env!("CARGO_PKG_VERSION_MINOR")) * 100 + p(env!("CARGO_PKG_VERSION_PATCH"))
}

#[no_mangle]
pub extern "C" fn getVersionStr() -> *mut c_char {
    text(VERSION)
}

#[no_mangle]
pub extern "C" fn getVersionEx() -> *mut c_char {
    text(&version_ex())
}

#[no_mangle]
pub extern "C" fn getExtendedAPIInfo() -> *mut c_char {
    text(&format!("roadrunner_nitro C API {}\nC API version: 1.0\nlibSBML version: none\n", version_ex()))
}

#[no_mangle]
pub extern "C" fn getBuildDate() -> *mut c_char {
    text(env!("RRC_BUILD_DATE"))
}

#[no_mangle]
pub extern "C" fn getBuildTime() -> *mut c_char {
    text(env!("RRC_BUILD_TIME"))
}

#[no_mangle]
pub extern "C" fn getBuildDateTime() -> *mut c_char {
    text(concat!(env!("RRC_BUILD_DATE"), " ", env!("RRC_BUILD_TIME")))
}

#[no_mangle]
pub extern "C" fn getCopyright() -> *mut c_char {
    text("roadrunner_nitro (tellurium-nitro). C API headers from libroadrunner, (C) University of Washington, Apache License 2.0")
}

#[no_mangle]
pub extern "C" fn getlibSBMLVersion(_handle: RRHandle) -> *mut c_char {
    text("none")
}

// ---------------------------------------------------------------- folders & compiler (per instance, stored only)

static INSTALL_FOLDER: Mutex<String> = Mutex::new(String::new());

/// Process-wide, as upstream's gInstallFolder.
#[no_mangle]
pub unsafe extern "C" fn setInstallFolder(folder: *const c_char) -> bool {
    run(false, || { *INSTALL_FOLDER.lock().unwrap() = st(folder)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn setTempFolder(handle: RRHandle, folder: *const c_char) -> bool {
    with(handle, false, |r| { r.temp = st(folder)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getTempFolder(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| Ok(text(&r.temp)))
}

#[no_mangle]
pub extern "C" fn getWorkingDirectory() -> *mut c_char {
    match std::env::current_dir() {
        Ok(d) => text(&d.to_string_lossy()),
        Err(e) => { set_error(e.to_string()); null_mut() }
    }
}

/// Directory containing this shared library.
#[no_mangle]
pub extern "C" fn getRRCAPILocation() -> *mut c_char {
    let mut info: libc::Dl_info = unsafe { std::mem::zeroed() };
    if unsafe { libc::dladdr(getRRCAPILocation as *const libc::c_void, &mut info) } == 0 || info.dli_fname.is_null() {
        set_error("getRRCAPILocation: dladdr failed");
        return null_mut();
    }
    let path = unsafe { cstr(info.dli_fname) }.unwrap_or_default();
    let dir = std::path::Path::new(&path).parent().map(|p| p.to_string_lossy().into_owned()).unwrap_or_default();
    text(&dir)
}

// The Nitro engine compiles in process (Cranelift): no external compiler. The settings are kept per
// instance and reported back, as libroadrunner's LLVM build accepts and ignores them.
#[no_mangle]
pub unsafe extern "C" fn setCompiler(handle: RRHandle, fNameWithPath: *const c_char) -> bool {
    with(handle, false, |r| { r.compiler = st(fNameWithPath)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getCompiler(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| Ok(text(&r.compiler)))
}

#[no_mangle]
pub unsafe extern "C" fn setCompilerLocation(handle: RRHandle, folder: *const c_char) -> bool {
    with(handle, false, |r| { r.compiler_location = st(folder)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getCompilerLocation(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| Ok(text(&r.compiler_location)))
}

#[no_mangle]
pub unsafe extern "C" fn setSupportCodeFolder(handle: RRHandle, folder: *const c_char) -> bool {
    with(handle, false, |r| { r.support = st(folder)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getSupportCodeFolder(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| Ok(text(&r.support)))
}

/// Upstream ignores the mode and returns false; so do we.
#[no_mangle]
pub extern "C" fn setCodeGenerationMode(_handle: RRHandle, _mode: c_int) -> bool {
    false
}

// ---------------------------------------------------------------- errors

#[no_mangle]
pub extern "C" fn hasError() -> bool {
    has_error()
}

/// Thread-local; do not free. "No Error" when nothing was recorded.
#[no_mangle]
pub extern "C" fn getLastError() -> *mut c_char {
    last_error_ptr()
}

// ---------------------------------------------------------------- config keys (process-wide map)
// ponytail: no upstream defaults; unknown keys fail with -1 until the core registers its keys.

static CONFIG: Mutex<BTreeMap<String, f64>> = Mutex::new(BTreeMap::new());

unsafe fn config_set(key: *const c_char, v: f64) -> c_int {
    let Some(k) = cstr(key) else { return not_set("NULL config key") };
    CONFIG.lock().unwrap().insert(k, v);
    1
}

unsafe fn config_get(key: *const c_char) -> Option<f64> {
    let k = cstr(key)?;
    let v = CONFIG.lock().unwrap().get(&k).copied();
    if v.is_none() {
        set_error(format!("RoadRunner exception: unknown config key '{k}'"));
    }
    v
}

fn not_set(msg: &str) -> c_int {
    set_error(msg);
    -1
}

#[no_mangle]
pub unsafe extern "C" fn setConfigBool(key: *const c_char, value: c_int) -> c_int {
    config_set(key, (value != 0) as c_int as f64)
}

#[no_mangle]
pub unsafe extern "C" fn getConfigBool(key: *const c_char) -> c_int {
    config_get(key).map_or(-1, |v| (v != 0.0) as c_int)
}

#[no_mangle]
pub unsafe extern "C" fn setConfigInt(key: *const c_char, value: c_int) -> c_int {
    config_set(key, value as f64)
}

#[no_mangle]
pub unsafe extern "C" fn getConfigInt(key: *const c_char) -> c_int {
    config_get(key).map_or(-1, |v| v as c_int)
}

#[no_mangle]
pub unsafe extern "C" fn setConfigDouble(key: *const c_char, value: f64) -> c_int {
    config_set(key, value)
}

#[no_mangle]
pub unsafe extern "C" fn getConfigDouble(key: *const c_char) -> f64 {
    config_get(key).unwrap_or(-1.0)
}

/// NULL when no keys are set (as upstream for an empty list).
#[no_mangle]
pub extern "C" fn getListOfConfigKeys() -> *mut RRStringArray {
    let keys: Vec<String> = CONFIG.lock().unwrap().keys().cloned().collect();
    if keys.is_empty() { null_mut() } else { string_array(&keys) }
}

/// Setting::TypeId code -> name (static string, do not free).
#[no_mangle]
pub extern "C" fn solverTypeToString(code: c_int) -> *const c_char {
    const NAMES: [&std::ffi::CStr; 11] =
        [c"empty", c"string", c"bool", c"int32", c"uint32", c"int64", c"uint64", c"float", c"double", c"char", c"uchar"];
    NAMES.get(code as usize).map_or(c"<invalid>".as_ptr(), |s| s.as_ptr())
}

// ---------------------------------------------------------------- instance plumbing

/// What an RRHandle points to: the engine plus the C-side per-instance state.
pub struct Inst {
    pub rr: RoadRunner,
    temp: String,
    compiler: String,
    compiler_location: String,
    support: String,
    /// integrators instantiated so far (getNumInstantiatedIntegrators)
    used: Vec<&'static str>,
    /// getSimulationResultAsDoubleMatrix's matrix (owned here, valid until the next call or freeRRInstance)
    last: *mut RRDoubleMatrix,
}

impl Inst {
    fn new(temp: String, compiler: String) -> Inst {
        Inst { rr: RoadRunner::new(), temp, compiler, compiler_location: String::new(), support: String::new(), used: vec!["cvode"], last: null_mut() }
    }
    fn set_integrator(&mut self, name: &str) -> R<()> {
        self.rr.set_integrator(name)?;
        let n = self.rr.integrator().name;
        if !self.used.contains(&n) { self.used.push(n); }
        Ok(())
    }
}

impl Drop for Inst {
    fn drop(&mut self) { unsafe { crate::utilities::freeMatrix(self.last); } }
}

type R<T> = Result<T, String>;

const NO_MODEL: &str = "A model needs to be loaded before one can use this method";

fn panic_text(p: Box<dyn std::any::Any + Send>) -> String {
    p.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| p.downcast_ref::<String>().cloned()).unwrap_or_else(|| "panic".into())
}

/// Run `f` (no panic crosses the FFI); an error is recorded as upstream's catch macros do
/// ("RoadRunner exception: <what>\n") and `fail` returned.
fn run<T>(fail: T, f: impl FnOnce() -> R<T>) -> T {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(v)) => v,
        Ok(Err(e)) => { set_error(format!("RoadRunner exception: {e}\n")); fail }
        Err(p) => { set_error(format!("RoadRunner exception: internal error: {}\n", panic_text(p))); fail }
    }
}

unsafe fn inst<'a>(h: RRHandle) -> R<&'a mut Inst> {
    (h as *mut Inst).as_mut().ok_or_else(|| "Failed to cast to a valid RoadRunner handle".to_string())
}

unsafe fn with<T>(h: RRHandle, fail: T, f: impl FnOnce(&mut Inst) -> R<T>) -> T {
    run(fail, || f(inst(h)?))
}

/// As `with`, for calls that need a loaded model.
unsafe fn wm<T>(h: RRHandle, fail: T, f: impl FnOnce(&mut RoadRunner) -> R<T>) -> T {
    run(fail, || {
        let r = inst(h)?;
        if !r.rr.is_model_loaded() { return Err(NO_MODEL.into()); }
        f(&mut r.rr)
    })
}

/// A C string argument.
unsafe fn st(p: *const c_char) -> R<String> { cstr(p).ok_or_else(|| "NULL string argument".to_string()) }

/// Write an out-argument.
unsafe fn out<T>(p: *mut T, v: T) -> R<bool> {
    if p.is_null() { return Err("NULL output argument".into()); }
    p.write(v);
    Ok(true)
}

/// RRVector argument as a Vec.
unsafe fn vec_arg(v: *const RRVector) -> R<Vec<f64>> {
    let v = v.as_ref().ok_or("NULL vector argument")?;
    Ok(if v.Count <= 0 || v.Data.is_null() { vec![] } else { std::slice::from_raw_parts(v.Data, v.Count as usize).to_vec() })
}

/// `char**` argument as strings.
unsafe fn strs(p: *const *const c_char, n: c_int) -> R<Vec<String>> {
    if n <= 0 { return Ok(vec![]); }
    if p.is_null() { return Err("NULL string list argument".into()); }
    (0..n as usize).map(|i| st(*p.add(i))).collect()
}

/// Upstream's createList(StringList) when the list may be empty: NULL for no items.
fn list_or_null(v: &[String]) -> *mut RRStringArray { if v.is_empty() { null_mut() } else { string_array(v) } }

/// createMatrix: NULL (no error) for an empty matrix, as upstream.
fn mat(m: Matrix) -> *mut RRDoubleMatrix { if m.data.is_empty() { null_mut() } else { matrix(m.nr, m.nc, &m.data) } }

fn cdata_of(m: &Matrix, weights: Option<&[f64]>) -> *mut RRCData { cdata(m.nr, &m.data, weights, &m.cols) }

/// Selection list text, split at commas and spaces outside parentheses (upstream's StringList(list,
/// ", ") also splits "cc(S1, k1)" apart, which then fails to parse).
fn split_selections(s: &str) -> Vec<String> {
    let (mut v, mut cur, mut depth) = (vec![], String::new(), 0i32);
    for ch in s.chars() {
        match ch {
            '(' => { depth += 1; cur.push(ch) }
            ')' => { depth -= 1; cur.push(ch) }
            ',' | ' ' | '\t' if depth <= 0 => { if !cur.is_empty() { v.push(std::mem::take(&mut cur)); } }
            _ => cur.push(ch),
        }
    }
    if !cur.is_empty() { v.push(cur); }
    v
}

/// Nested RRList items (upstream's ArrayList).
enum L { S(String), L(Vec<L>) }

/// createArrayList: NULL for an empty list.
fn rrlist(items: &[L]) -> *mut RRList {
    if items.is_empty() { return null_mut(); }
    let l = crate::utilities::createRRList();
    for it in items {
        let mut p = match it {
            L::S(s) => { let c = std::ffi::CString::new(s.replace('\0', " ")).unwrap(); unsafe { crate::utilities::createStringItem(c.as_ptr()) } }
            L::L(v) => crate::utilities::createListItem(rrlist(v)),
        };
        unsafe { crate::utilities::addItem(l, &mut p) };
    }
    l
}

/// {label, {items}}
fn group(label: &str, items: Vec<String>) -> L { L::L(vec![L::S(label.into()), L::L(items.into_iter().map(L::S).collect())]) }

/// {{name, {items}}, ...}
fn groups(v: Vec<(String, Vec<String>)>) -> Vec<L> { v.into_iter().map(|(n, i)| group(&n, i)).collect() }

fn idx(i: c_int, n: usize, what: &str) -> R<usize> {
    usize::try_from(i).ok().filter(|&i| i < n).ok_or_else(|| format!("Index in {what} out of range: [{i}]"))
}

fn conc(ids: Vec<String>) -> Vec<String> { ids.into_iter().map(|s| format!("[{s}]")).collect() }

// ---------------------------------------------------------------- instances

#[no_mangle]
pub extern "C" fn createRRInstance() -> RRHandle {
    run(null_mut(), || Ok(Box::into_raw(Box::new(Inst::new(std::env::temp_dir().to_string_lossy().into_owned(), String::new()))) as RRHandle))
}

#[no_mangle]
pub unsafe extern "C" fn createRRInstanceEx(tempFolder: *const c_char, compiler: *const c_char) -> RRHandle {
    run(null_mut(), || {
        let temp = match cstr(tempFolder) {
            Some(t) if !std::path::Path::new(&t).exists() => return Err(format!("The temporary folder: {t} do not exist")),
            Some(t) => t,
            None => std::env::temp_dir().to_string_lossy().into_owned(),
        };
        Ok(Box::into_raw(Box::new(Inst::new(temp, cstr(compiler).unwrap_or_default()))) as RRHandle)
    })
}

#[no_mangle]
pub unsafe extern "C" fn freeRRInstance(handle: RRHandle) -> bool {
    run(false, || { inst(handle)?; drop(Box::from_raw(handle as *mut Inst)); Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getInfo(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| Ok(text(&format!(
        "roadrunner_nitro C API {}\nModel loaded: {}\nIntegrator: {}\nSteady state solver: {}\nTemp folder: {}\n",
        version_ex(), r.rr.is_model_loaded(), r.rr.integrator().name, r.rr.steady_state_solver().name, r.temp))))
}

// ---------------------------------------------------------------- flags

#[no_mangle]
pub unsafe extern "C" fn setComputeAndAssignConservationLaws(handle: RRHandle, On_Or_Off: bool) -> bool {
    with(handle, false, |r| { r.rr.set_conservation_laws(On_Or_Off)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn getComputeAndAssignConservationLaws(handle: RRHandle, value: *mut c_int) -> c_int {
    with(handle, 0, |r| out(value, r.rr.conservation_laws() as c_int).map(c_int::from))
}

// ---------------------------------------------------------------- load / save

unsafe fn load_file(handle: RRHandle, fileName: *const c_char) -> bool {
    with(handle, false, |r| {
        let f = st(fileName)?;
        if !std::path::Path::new(&f).exists() {
            set_error(format!("The file {f} was not found")); // upstream: no "RoadRunner exception" prefix
            return Ok(false);
        }
        r.rr.load_sbml(&f, None)?;
        Ok(true)
    })
}

#[no_mangle]
pub unsafe extern "C" fn loadSBMLFromFile(handle: RRHandle, fileName: *const c_char) -> bool { load_file(handle, fileName) }

/// forceRecompile: the engine compiles every load.
#[no_mangle]
pub unsafe extern "C" fn loadSBMLFromFileE(handle: RRHandle, fileName: *const c_char, _forceRecompile: bool) -> bool { load_file(handle, fileName) }

#[no_mangle]
pub unsafe extern "C" fn loadSBML(handle: RRHandle, sbml: *const c_char) -> bool {
    with(handle, false, |r| { r.rr.load_sbml(&st(sbml)?, None)?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn loadSBMLEx(handle: RRHandle, sbml: *const c_char, _forceRecompile: bool) -> bool { loadSBML(handle, sbml) }

#[no_mangle]
pub unsafe extern "C" fn saveState(handle: RRHandle, filename: *const c_char) -> bool {
    with(handle, false, |r| { r.rr.save_state(std::path::Path::new(&st(filename)?))?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn loadState(handle: RRHandle, filename: *const c_char) -> bool {
    with(handle, false, |r| { r.rr.load_state(std::path::Path::new(&st(filename)?))?; Ok(true) })
}

/// True if a model was loaded (and is now cleared).
#[no_mangle]
pub unsafe extern "C" fn clearModel(handle: RRHandle) -> bool {
    with(handle, false, |r| { let was = r.rr.is_model_loaded(); r.rr.clear_model(); Ok(was) })
}

#[no_mangle]
pub unsafe extern "C" fn isModelLoaded(handle: RRHandle) -> bool {
    with(handle, false, |r| Ok(r.rr.is_model_loaded()))
}

/// SBML Test Suite settings file: start, duration, steps -> the time course; absolute, relative ->
/// the current integrator's tolerances (upstream applies the time settings to a copy, i.e. not at all).
#[no_mangle]
pub unsafe extern "C" fn loadSimulationSettings(handle: RRHandle, fileName: *const c_char) -> bool {
    with(handle, false, |r| {
        let f = st(fileName)?;
        if !std::path::Path::new(&f).exists() {
            set_error(format!("The file {f} was not found"));
            return Ok(false);
        }
        let text = std::fs::read_to_string(&f).map_err(|e| format!("{f}: {e}"))?;
        let kv: BTreeMap<String, String> = text.lines().filter_map(|l| l.split_once(':')).map(|(k, v)| (k.trim().to_string(), v.trim().to_string())).collect();
        let num = |k: &str| -> R<Option<f64>> {
            kv.get(k).filter(|v| !v.is_empty()).map(|v| v.parse::<f64>().map_err(|_| format!("{f}: bad {k} '{v}'"))).transpose()
        };
        let rr = &mut r.rr;
        if let Some(s) = num("start")? { rr.end += s - rr.start; rr.start = s; rr.times = None; }
        if let Some(d) = num("duration")? { rr.end = rr.start + d; rr.times = None; }
        if let Some(n) = num("steps")? { rr.points = n as usize + 1; rr.times = None; }
        for (k, key) in [("absolute", "absolute_tolerance"), ("relative", "relative_tolerance")] {
            if let (Some(v), true) = (num(k)?, rr.integrator().get(key).is_ok()) { rr.integrator_mut().set(key, Value::Double(v))?; }
        }
        Ok(true)
    })
}

#[no_mangle]
pub unsafe extern "C" fn getCurrentSBML(handle: RRHandle) -> *mut c_char {
    wm(handle, null_mut(), |rr| Ok(text(&rr.get_current_sbml()?)))
}

#[no_mangle]
pub unsafe extern "C" fn getSBML(handle: RRHandle) -> *mut c_char {
    wm(handle, null_mut(), |rr| Ok(text(&rr.get_sbml()?)))
}

// ---------------------------------------------------------------- model editing

/// Apply an SBML text edit; `regen`: recompile (the NoRegen variants defer to regenerate()).
unsafe fn ed(handle: RRHandle, regen: bool, f: impl FnOnce(&str) -> R<String>) -> bool {
    with(handle, false, |r| { r.rr.edit(|s| { *s = f(s)?; Ok(()) }, regen)?; Ok(true) })
}

macro_rules! edits {
    ($($name:ident, $noregen:ident ($($arg:ident: $t:ty),*) => |$s:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle, $($arg: $t),*) -> bool { ed(handle, true, |$s| $body) }
        #[no_mangle]
        pub unsafe extern "C" fn $noregen(handle: RRHandle, $($arg: $t),*) -> bool { ed(handle, false, |$s| $body) }
    )*};
}

edits! {
    addSpeciesAmount, addSpeciesAmountNoRegen(sid: *const c_char, compartment: *const c_char, initialAmount: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool)
        => |s| edit::add_species_amount(s, &st(sid)?, &st(compartment)?, initialAmount, hasOnlySubstanceUnits, boundaryCondition, "");
    addSpeciesConcentration, addSpeciesConcentrationNoRegen(sid: *const c_char, compartment: *const c_char, initialConcentration: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool)
        => |s| edit::add_species_concentration(s, &st(sid)?, &st(compartment)?, initialConcentration, hasOnlySubstanceUnits, boundaryCondition, "");
    removeSpecies, removeSpeciesNoRegen(sid: *const c_char) => |s| edit::remove_species(s, &st(sid)?);
    setBoundary, setBoundaryNoRegen(sid: *const c_char, boundaryCondition: bool) => |s| edit::set_boundary(s, &st(sid)?, boundaryCondition);
    setHasOnlySubstanceUnits, setHasOnlySubstanceUnitsNoRegen(sid: *const c_char, hasOnlySubstanceUnits: bool)
        => |s| edit::set_has_only_substance_units(s, &st(sid)?, hasOnlySubstanceUnits);
    setConstant, setConstantNoRegen(sid: *const c_char, constant: bool) => |s| edit::set_constant(s, &st(sid)?, constant);
    addReactionFromSBML, addReactionFromSBMLNoRegen(sbmlRep: *const c_char) => |s| edit::add_reaction_from_sbml(s, &st(sbmlRep)?);
    removeReaction, removeReactionNoRegen(rid: *const c_char) => |s| edit::remove_reaction(s, &st(rid)?, false);
    setReversible, setReversibleNoRegen(rid: *const c_char, reversible: bool) => |s| edit::set_reversible(s, &st(rid)?, reversible);
    setKineticLaw, setKineticLawNoRegen(rid: *const c_char, kineticLaw: *const c_char) => |s| edit::set_kinetic_law(s, &st(rid)?, &st(kineticLaw)?);
    addParameter, addParameterNoRegen(pid: *const c_char, value: f64) => |s| edit::add_parameter(s, &st(pid)?, value);
    removeParameter, removeParameterNoRegen(pid: *const c_char) => |s| edit::remove_parameter(s, &st(pid)?);
    addCompartment, addCompartmentNoRegen(cid: *const c_char, initVolume: f64) => |s| edit::add_compartment(s, &st(cid)?, initVolume);
    removeCompartment, removeCompartmentNoRegen(cid: *const c_char) => |s| edit::remove_compartment(s, &st(cid)?);
    addAssignmentRule, addAssignmentRuleNoRegen(vid: *const c_char, formula: *const c_char) => |s| edit::add_assignment_rule(s, &st(vid)?, &st(formula)?);
    addRateRule, addRateRuleNoRegen(vid: *const c_char, formula: *const c_char) => |s| edit::add_rate_rule(s, &st(vid)?, &st(formula)?);
    removeRules, removeRulesNoRegen(vid: *const c_char) => |s| edit::remove_rules(s, &st(vid)?);
    addEvent, addEventNoRegen(eid: *const c_char, useValuesFromTriggerTime: bool, trigger: *const c_char)
        => |s| edit::add_event(s, &st(eid)?, useValuesFromTriggerTime, &st(trigger)?);
    addTrigger, addTriggerNoRegen(eid: *const c_char, trigger: *const c_char) => |s| edit::add_trigger(s, &st(eid)?, &st(trigger)?);
    setPersistent, setPersistentNoRegen(eid: *const c_char, persistent: bool) => |s| edit::set_persistent(s, &st(eid)?, persistent);
    setTriggerInitialValue, setTriggerInitialValueNoRegen(eid: *const c_char, initValue: bool) => |s| edit::set_trigger_initial_value(s, &st(eid)?, initValue);
    addPriority, addPriorityNoRegen(eid: *const c_char, priority: *const c_char) => |s| edit::add_priority(s, &st(eid)?, &st(priority)?);
    addDelay, addDelayNoRegen(eid: *const c_char, delay: *const c_char) => |s| edit::add_delay(s, &st(eid)?, &st(delay)?);
    addEventAssignment, addEventAssignmentNoRegen(eid: *const c_char, vid: *const c_char, formula: *const c_char)
        => |s| edit::add_event_assignment(s, &st(eid)?, &st(vid)?, &st(formula)?);
    removeEventAssignments, removeEventAssignmentsNoRegen(eid: *const c_char, vid: *const c_char) => |s| edit::remove_event_assignments(s, &st(eid)?, &st(vid)?);
    removeEvent, removeEventsNoRegen(eid: *const c_char) => |s| edit::remove_event(s, &st(eid)?);
}

#[no_mangle]
pub unsafe extern "C" fn setInitAmount(handle: RRHandle, sid: *const c_char, initAmount: f64) -> bool {
    ed(handle, true, |s| edit::set_init_amount(s, &st(sid)?, initAmount))
}

#[no_mangle]
pub unsafe extern "C" fn setInitConcentration(handle: RRHandle, sid: *const c_char, initConcentration: f64) -> bool {
    ed(handle, true, |s| edit::set_init_concentration(s, &st(sid)?, initConcentration))
}

#[no_mangle]
pub unsafe extern "C" fn addReaction(handle: RRHandle, rid: *const c_char, reactants: *mut *const c_char, numReactants: c_int,
                                     products: *mut *const c_char, numProducts: c_int, kineticLaw: *const c_char) -> bool {
    add_reaction(handle, true, rid, reactants, numReactants, products, numProducts, kineticLaw)
}

#[no_mangle]
pub unsafe extern "C" fn addReactionNoRegen(handle: RRHandle, rid: *const c_char, reactants: *mut *const c_char, numReactants: c_int,
                                            products: *mut *const c_char, numProducts: c_int, kineticLaw: *const c_char) -> bool {
    add_reaction(handle, false, rid, reactants, numReactants, products, numProducts, kineticLaw)
}

#[allow(clippy::too_many_arguments)]
unsafe fn add_reaction(handle: RRHandle, regen: bool, rid: *const c_char, re: *mut *const c_char, nre: c_int,
                       pr: *mut *const c_char, npr: c_int, law: *const c_char) -> bool {
    ed(handle, regen, |s| {
        let (re, pr) = (strs(re, nre)?, strs(pr, npr)?);
        let re: Vec<&str> = re.iter().map(String::as_str).collect();
        let pr: Vec<&str> = pr.iter().map(String::as_str).collect();
        edit::add_reaction(s, &st(rid)?, &re, &pr, &st(law)?)
    })
}

#[no_mangle]
pub unsafe extern "C" fn addInitialAssignment(handle: RRHandle, vid: *const c_char, formula: *const c_char, forceRegenerate: bool) -> bool {
    ed(handle, forceRegenerate, |s| edit::add_initial_assignment(s, &st(vid)?, &st(formula)?))
}

#[no_mangle]
pub unsafe extern "C" fn getKineticLaw(handle: RRHandle, rid: *const c_char) -> *mut c_char {
    wm(handle, null_mut(), |rr| Ok(text(&edit::get_kinetic_law(&rr.get_sbml()?, &st(rid)?)?)))
}

#[no_mangle]
pub unsafe extern "C" fn regenerate(handle: RRHandle) -> bool {
    with(handle, false, |r| { r.rr.regenerate()?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn validateCurrentSBML(handle: RRHandle) -> bool {
    wm(handle, false, |rr| { edit::validate_current_sbml(&rr.get_sbml()?)?; Ok(true) })
}

// ---------------------------------------------------------------- SBML utilities, configuration XML

/// `sArg`: SBML text, or the path of an SBML file (as upstream).
#[no_mangle]
pub unsafe extern "C" fn getParamPromotedSBML(handle: RRHandle, sArg: *const c_char) -> *mut c_char {
    with(handle, null_mut(), |_| {
        let a = st(sArg)?;
        let sbml = if a.trim_start().starts_with('<') { a } else { std::fs::read_to_string(&a).map_err(|e| format!("{a}: {e}"))? };
        Ok(text(&edit::get_param_promoted_sbml(&sbml)?))
    })
}

fn type_code(v: &Value) -> c_int {
    match v { Value::Str(_) => 1, Value::Bool(_) => 2, Value::Int(_) => 3, Value::UInt(_) => 6, Value::Double(_) => 8, Value::Doubles(_) => 11 }
}

/// Setting::getValueAsString.
fn value_text(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        Value::Bool(b) => b.to_string(),
        Value::Int(i) => i.to_string(),
        Value::UInt(u) => u.to_string(),
        Value::Double(d) => crate::utilities::g(*d),
        Value::Doubles(d) => format!("[{}]", d.iter().map(|x| crate::utilities::g(*x)).collect::<Vec<_>>().join(", ")),
    }
}

fn xml_escape(s: &str) -> String { s.replace('&', "&amp;").replace('"', "&quot;").replace('<', "&lt;").replace('>', "&gt;") }

/// The current integrator's and steady-state solver's settings in the caps format of rrc_api.h.
#[no_mangle]
pub unsafe extern "C" fn getConfigurationXML(handle: RRHandle) -> *mut c_char {
    with(handle, null_mut(), |r| {
        let mut x = String::from("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<caps name=\"RoadRunner\" description=\"Settings For RoadRunner\">\n");
        for (sec, s) in [("Integration", r.rr.integrator()), ("SteadyState", r.rr.steady_state_solver())] {
            x += &format!("  <section name=\"{sec}\" method=\"{}\" description=\"{}\">\n", s.name, xml_escape(s.description));
            for p in &s.params {
                let ty = match &p.value { Value::Str(_) => "string", Value::Bool(_) => "bool", Value::Int(_) | Value::UInt(_) => "integer", Value::Double(_) => "double", Value::Doubles(_) => "doubles" };
                x += &format!("    <cap name=\"{}\" value=\"{}\" hint=\"{}\" type=\"{ty}\" />\n", p.name, xml_escape(&value_text(&p.value)), xml_escape(p.hint));
            }
            x += "  </section>\n";
        }
        x += "</caps>\n";
        Ok(text(&x))
    })
}

fn xml_attr(tag: &str, name: &str) -> Option<String> {
    let k = format!("{name}=\"");
    let i = tag.find(&k)? + k.len();
    let v = &tag[i..i + tag[i..].find('"')?];
    Some(v.replace("&quot;", "\"").replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&"))
}

/// Apply a caps document (getConfigurationXML's format): each section's method selects the
/// integrator / steady-state solver, each cap sets one of its settings.
#[no_mangle]
pub unsafe extern "C" fn setConfigurationXML(handle: RRHandle, caps: *const c_char) -> bool {
    with(handle, false, |r| {
        let xml = st(caps)?;
        let mut ss = false;
        for tag in xml.split('<').skip(1).map(|t| t.split('>').next().unwrap_or("")) {
            if tag.starts_with("section") {
                ss = xml_attr(tag, "name").is_some_and(|n| n.eq_ignore_ascii_case("SteadyState"));
                if let Some(m) = xml_attr(tag, "method") {
                    let m = m.to_lowercase();
                    if ss { r.rr.set_steady_state_solver(&m)? } else { r.set_integrator(&m)? }
                }
            } else if tag.starts_with("cap ") {
                let (Some(n), Some(v)) = (xml_attr(tag, "name"), xml_attr(tag, "value")) else { return Err(format!("cap without name or value: <{tag}>")) };
                let val = match v.as_str() {
                    "true" => Value::Bool(true),
                    "false" => Value::Bool(false),
                    _ if v.starts_with('[') => Value::Doubles(v.trim_matches(['[', ']']).split(',').filter(|x| !x.trim().is_empty())
                        .map(|x| x.trim().parse::<f64>().map_err(|_| format!("cap {n}: bad value '{v}'"))).collect::<R<_>>()?),
                    _ => v.parse::<f64>().map(Value::Double).unwrap_or(Value::Str(v.clone())),
                };
                let s = if ss { r.rr.steady_state_solver_mut() } else { r.rr.integrator_mut() };
                s.set(&n, val)?;
            }
        }
        Ok(true)
    })
}

// ---------------------------------------------------------------- integrators & steady-state solvers

fn registered(names: &[&'static str], n: c_int, ss: bool, f: impl Fn(&Solver) -> String) -> *mut c_char {
    run(null_mut(), || {
        let name = *names.get(n.max(0) as usize).ok_or_else(|| format!("No solver registered at index {n}"))?;
        let s = if ss { Solver::steady_state(name) } else { Solver::integrator(name) }.ok_or("unknown solver")?;
        Ok(text(&f(&s)))
    })
}

#[no_mangle]
pub extern "C" fn getNumRegisteredIntegrators() -> c_int { INTEGRATORS.len() as c_int }

#[no_mangle]
pub extern "C" fn getRegisteredIntegratorName(n: c_int) -> *mut c_char { registered(&INTEGRATORS, n, false, |s| s.name.into()) }

#[no_mangle]
pub extern "C" fn getRegisteredIntegratorHint(n: c_int) -> *mut c_char { registered(&INTEGRATORS, n, false, |s| s.hint.into()) }

#[no_mangle]
pub extern "C" fn getRegisteredIntegratorDescription(n: c_int) -> *mut c_char { registered(&INTEGRATORS, n, false, |s| s.description.into()) }

#[no_mangle]
pub extern "C" fn getNumRegisteredSteadyStateSolvers() -> c_int { STEADY_STATE_SOLVERS.len() as c_int }

#[no_mangle]
pub extern "C" fn getRegisteredSteadyStateSolverName(n: c_int) -> *mut c_char { registered(&STEADY_STATE_SOLVERS, n, true, |s| s.name.into()) }

#[no_mangle]
pub extern "C" fn getRegisteredSteadyStateSolverHint(n: c_int) -> *mut c_char { registered(&STEADY_STATE_SOLVERS, n, true, |s| s.hint.into()) }

#[no_mangle]
pub extern "C" fn getRegisteredSteadyStateSolverDescription(n: c_int) -> *mut c_char { registered(&STEADY_STATE_SOLVERS, n, true, |s| s.description.into()) }

#[no_mangle]
pub unsafe extern "C" fn getNumInstantiatedIntegrators(handle: RRHandle) -> c_int {
    with(handle, -1, |r| Ok(r.used.len() as c_int))
}

#[no_mangle]
pub unsafe extern "C" fn setCurrentIntegrator(handle: RRHandle, nameOfIntegrator: *const c_char) -> c_int {
    with(handle, -1, |r| { r.set_integrator(&st(nameOfIntegrator)?)?; Ok(1) })
}

#[no_mangle]
pub unsafe extern "C" fn setCurrentSteadyStateSolver(handle: RRHandle, nameOfSteadyStateSolver: *const c_char) -> c_int {
    with(handle, -1, |r| { r.rr.set_steady_state_solver(&st(nameOfSteadyStateSolver)?)?; Ok(1) })
}

/// CVODEIntegrator::setIndividualTolerance: absolute tolerance on the amount of one floating species
/// or rate-rule variable; absolute_tolerance becomes a vector (getCurrentIntegratorParameterDoubleArray).
#[no_mangle]
pub unsafe extern "C" fn setCurrentIntegratorIndividualTolerance(handle: RRHandle, sid: *const c_char, value: f64) -> c_int {
    wm(handle, -1, |rr| { rr.set_individual_tolerance(&st(sid)?, value)?; Ok(1) })
}

fn integ(rr: &mut RoadRunner) -> &mut Solver { rr.integrator_mut() }
fn ssolv(rr: &mut RoadRunner) -> &mut Solver { rr.steady_state_solver_mut() }

/// (name, display name, hint, value) of setting `n`.
type P<'a> = (&'static str, &'static str, &'static str, &'a Value);

fn nth(s: &Solver, n: c_int) -> R<P<'_>> {
    let p = usize::try_from(n).ok().and_then(|i| s.params.get(i)).ok_or_else(|| format!("{}: no setting at index {n}", s.name))?;
    Ok((p.name, p.display, p.hint, &p.value))
}

fn param<'a>(s: &'a Solver, name: &str) -> R<P<'a>> {
    let p = s.params.iter().find(|p| p.name == name).ok_or_else(|| format!("{}: no setting '{name}'", s.name))?;
    Ok((p.name, p.display, p.hint, &p.value))
}

fn num(v: &Value, what: &str) -> R<f64> {
    match v { Value::Str(_) | Value::Doubles(_) => Err(format!("setting '{what}' is not a number")), v => Ok(v.f64()) }
}

macro_rules! solver_api {
    ($sel:ident; $name:ident $desc:ident $hint:ident $count:ident $nname:ident $ndisp:ident $ndesc:ident $nhint:ident $ntype:ident
     $reset:ident $names:ident $pdesc:ident $phint:ident $ptype:ident $gi:ident $si:ident $gu:ident $su:ident $gd:ident $sd:ident
     $gs:ident $ss:ident $gb:ident $sb:ident) => {
        #[no_mangle] pub unsafe extern "C" fn $name(handle: RRHandle) -> *mut c_char { with(handle, null_mut(), |r| Ok(text($sel(&mut r.rr).name))) }
        #[no_mangle] pub unsafe extern "C" fn $desc(handle: RRHandle) -> *mut c_char { with(handle, null_mut(), |r| Ok(text($sel(&mut r.rr).description))) }
        #[no_mangle] pub unsafe extern "C" fn $hint(handle: RRHandle) -> *mut c_char { with(handle, null_mut(), |r| Ok(text($sel(&mut r.rr).hint))) }
        #[no_mangle] pub unsafe extern "C" fn $count(handle: RRHandle) -> c_int { with(handle, -1, |r| Ok($sel(&mut r.rr).params.len() as c_int)) }
        #[no_mangle] pub unsafe extern "C" fn $nname(handle: RRHandle, n: c_int) -> *mut c_char { with(handle, null_mut(), |r| Ok(text(nth($sel(&mut r.rr), n)?.0))) }
        #[no_mangle] pub unsafe extern "C" fn $ndisp(handle: RRHandle, n: c_int) -> *mut c_char { with(handle, null_mut(), |r| Ok(text(nth($sel(&mut r.rr), n)?.1))) }
        // ponytail: the core keeps one help text per setting (its hint); it serves as the description too
        #[no_mangle] pub unsafe extern "C" fn $ndesc(handle: RRHandle, n: c_int) -> *mut c_char { with(handle, null_mut(), |r| Ok(text(nth($sel(&mut r.rr), n)?.2))) }
        #[no_mangle] pub unsafe extern "C" fn $nhint(handle: RRHandle, n: c_int) -> *mut c_char { with(handle, null_mut(), |r| Ok(text(nth($sel(&mut r.rr), n)?.2))) }
        #[no_mangle] pub unsafe extern "C" fn $ntype(handle: RRHandle, n: c_int) -> c_int { with(handle, -1, |r| Ok(type_code(&nth($sel(&mut r.rr), n)?.3))) }
        #[no_mangle] pub unsafe extern "C" fn $reset(handle: RRHandle) -> c_int { with(handle, -1, |r| { $sel(&mut r.rr).reset(); Ok(1) }) }
        #[no_mangle] pub unsafe extern "C" fn $names(handle: RRHandle) -> *mut RRStringArray {
            with(handle, null_mut(), |r| Ok(string_array(&$sel(&mut r.rr).params.iter().map(|p| p.name).collect::<Vec<_>>())))
        }
        #[no_mangle] pub unsafe extern "C" fn $pdesc(handle: RRHandle, parameterName: *const c_char) -> *mut c_char {
            with(handle, null_mut(), |r| Ok(text(param($sel(&mut r.rr), &st(parameterName)?)?.2)))
        }
        #[no_mangle] pub unsafe extern "C" fn $phint(handle: RRHandle, parameterName: *const c_char) -> *mut c_char {
            with(handle, null_mut(), |r| Ok(text(param($sel(&mut r.rr), &st(parameterName)?)?.2)))
        }
        #[no_mangle] pub unsafe extern "C" fn $ptype(handle: RRHandle, parameterName: *const c_char) -> c_int {
            with(handle, -1, |r| Ok(type_code(&param($sel(&mut r.rr), &st(parameterName)?)?.3)))
        }
        #[no_mangle] pub unsafe extern "C" fn $gi(handle: RRHandle, parameterName: *const c_char) -> c_int {
            with(handle, -1, |r| { let k = st(parameterName)?; Ok(num($sel(&mut r.rr).get(&k)?, &k)? as c_int) })
        }
        #[no_mangle] pub unsafe extern "C" fn $si(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int {
            with(handle, -1, |r| { $sel(&mut r.rr).set(&st(parameterName)?, Value::Int(value as i64))?; Ok(1) })
        }
        #[no_mangle] pub unsafe extern "C" fn $gu(handle: RRHandle, parameterName: *const c_char) -> c_uint {
            with(handle, u32::MAX, |r| { let k = st(parameterName)?; Ok(num($sel(&mut r.rr).get(&k)?, &k)? as c_uint) })
        }
        #[no_mangle] pub unsafe extern "C" fn $su(handle: RRHandle, parameterName: *const c_char, value: c_uint) -> c_int {
            with(handle, -1, |r| { $sel(&mut r.rr).set(&st(parameterName)?, Value::UInt(value as u64))?; Ok(1) })
        }
        #[no_mangle] pub unsafe extern "C" fn $gd(handle: RRHandle, parameterName: *const c_char) -> f64 {
            with(handle, -1.0, |r| { let k = st(parameterName)?; num($sel(&mut r.rr).get(&k)?, &k) })
        }
        #[no_mangle] pub unsafe extern "C" fn $sd(handle: RRHandle, parameterName: *const c_char, value: f64) -> c_int {
            with(handle, -1, |r| { $sel(&mut r.rr).set(&st(parameterName)?, Value::Double(value))?; Ok(1) })
        }
        #[no_mangle] pub unsafe extern "C" fn $gs(handle: RRHandle, parameterName: *const c_char) -> *mut c_char {
            with(handle, null_mut(), |r| Ok(text(&value_text($sel(&mut r.rr).get(&st(parameterName)?)?))))
        }
        #[no_mangle] pub unsafe extern "C" fn $ss(handle: RRHandle, parameterName: *const c_char, value: *const c_char) -> c_int {
            with(handle, -1, |r| { $sel(&mut r.rr).set(&st(parameterName)?, Value::Str(st(value)?))?; Ok(1) })
        }
        #[no_mangle] pub unsafe extern "C" fn $gb(handle: RRHandle, parameterName: *const c_char) -> c_int {
            with(handle, -1, |r| { let k = st(parameterName)?; Ok((num($sel(&mut r.rr).get(&k)?, &k)? != 0.0) as c_int) })
        }
        #[no_mangle] pub unsafe extern "C" fn $sb(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int {
            with(handle, -1, |r| { $sel(&mut r.rr).set(&st(parameterName)?, Value::Bool(value != 0))?; Ok(1) })
        }
    };
}

solver_api!(integ; getCurrentIntegratorName getCurrentIntegratorDescription getCurrentIntegratorHint getNumberOfCurrentIntegratorParameters
    getCurrentIntegratorNthParameterName getCurrentIntegratorNthParameterDisplayName getCurrentIntegratorNthParameterDescription
    getCurrentIntegratorNthParameterHint getCurrentIntegratorNthParameterType resetCurrentIntegratorParameters
    getListOfCurrentIntegratorParameterNames getCurrentIntegratorParameterDescription getCurrentIntegratorParameterHint
    getCurrentIntegratorParameterType getCurrentIntegratorParameterInt setCurrentIntegratorParameterInt
    getCurrentIntegratorParameterUInt setCurrentIntegratorParameterUInt getCurrentIntegratorParameterDouble
    setCurrentIntegratorParameterDouble getCurrentIntegratorParameterString setCurrentIntegratorParameterString
    getCurrentIntegratorParameterBoolean setCurrentIntegratorParameterBoolean);

solver_api!(ssolv; getCurrentSteadyStateSolverName getCurrentSteadyStateSolverDescription getCurrentSteadyStateSolverHint
    getNumberOfCurrentSteadyStateSolverParameters getCurrentSteadyStateSolverNthParameterName
    getCurrentSteadyStateSolverNthParameterDisplayName getCurrentSteadyStateSolverNthParameterDescription
    getCurrentSteadyStateSolverNthParameterHint getCurrentSteadyStateSolverNthParameterType resetCurrentSteadyStateSolverParameters
    getListOfCurrentSteadyStateSolverParameterNames getCurrentSteadyStateSolverParameterDescription
    getCurrentSteadyStateSolverParameterHint getCurrentSteadyStateSolverParameterType getCurrentSteadyStateSolverParameterInt
    setCurrentSteadyStateSolverParameterInt getCurrentSteadyStateSolverParameterUInt setCurrentSteadyStateSolverParameterUInt
    getCurrentSteadyStateSolverParameterDouble setCurrentSteadyStateSolverParameterDouble getCurrentSteadyStateSolverParameterString
    setCurrentSteadyStateSolverParameterString getCurrentSteadyStateSolverParameterBoolean setCurrentSteadyStateSolverParameterBoolean);

/// `*value`: malloc'd copy (free with free()), as upstream.
#[no_mangle]
pub unsafe extern "C" fn getCurrentIntegratorParameterDoubleArray(handle: RRHandle, parameterName: *const c_char, value: *mut *mut f64, len: *mut c_int) -> c_int {
    with(handle, -1, |r| {
        let k = st(parameterName)?;
        let Value::Doubles(v) = r.rr.integrator().get(&k)? else { return Err(format!("setting '{k}' is not a double array")) };
        if value.is_null() || len.is_null() { return Err("NULL output argument".into()); }
        let p = libc::malloc(v.len().max(1) * std::mem::size_of::<f64>()) as *mut f64;
        if p.is_null() { return Err("out of memory".into()); }
        p.copy_from_nonoverlapping(v.as_ptr(), v.len());
        (*value, *len) = (p, v.len() as c_int);
        Ok(1)
    })
}

#[no_mangle]
pub unsafe extern "C" fn setCurrentIntegratorParameterDoubleArray(handle: RRHandle, parameterName: *const c_char, value: *mut f64, len: c_int) -> c_int {
    with(handle, -1, |r| {
        let v = if len <= 0 { vec![] } else if value.is_null() { return Err("NULL array argument".into()) } else { std::slice::from_raw_parts(value, len as usize).to_vec() };
        r.rr.integrator_mut().set(&st(parameterName)?, Value::Doubles(v))?;
        Ok(1)
    })
}

// ---------------------------------------------------------------- time course settings

/// Upstream keeps start and duration: a new start moves the end.
#[no_mangle]
pub unsafe extern "C" fn setTimeStart(handle: RRHandle, timeStart: f64) -> bool {
    with(handle, false, |r| { let rr = &mut r.rr; rr.end += timeStart - rr.start; rr.start = timeStart; rr.times = None; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn setTimeEnd(handle: RRHandle, timeEnd: f64) -> bool {
    with(handle, false, |r| { r.rr.end = timeEnd; r.rr.times = None; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn setNumPoints(handle: RRHandle, numberOfPoints: c_int) -> bool {
    with(handle, false, |r| {
        if numberOfPoints < 1 { return Err(format!("number of points must be at least 1, got {numberOfPoints}")); }
        r.rr.points = numberOfPoints as usize;
        r.rr.times = None;
        Ok(true)
    })
}

#[no_mangle]
pub unsafe extern "C" fn setTimes(handle: RRHandle, times: *const f64, size: c_int) -> bool {
    with(handle, false, |r| {
        if size <= 0 || times.is_null() { return Err("setTimes: no times".into()); }
        let t = std::slice::from_raw_parts(times, size as usize).to_vec();
        (r.rr.start, r.rr.end, r.rr.points) = (t[0], t[t.len() - 1], t.len());
        r.rr.times = Some(t);
        Ok(true)
    })
}

#[no_mangle]
pub unsafe extern "C" fn getTimeStart(handle: RRHandle, timeStart: *mut f64) -> bool { with(handle, false, |r| out(timeStart, r.rr.start)) }

#[no_mangle]
pub unsafe extern "C" fn getTimeEnd(handle: RRHandle, timeEnd: *mut f64) -> bool { with(handle, false, |r| out(timeEnd, r.rr.end)) }

#[no_mangle]
pub unsafe extern "C" fn getNumPoints(handle: RRHandle, numPoints: *mut c_int) -> bool { with(handle, false, |r| out(numPoints, r.rr.points as c_int)) }

unsafe fn set_selections(handle: RRHandle, steady: bool, sel: impl FnOnce() -> R<Vec<String>>) -> bool {
    wm(handle, false, |rr| {
        let s = sel()?;
        if steady { rr.set_steady_state_selections(&s)? } else { rr.set_selections(&s)? }
        Ok(true)
    })
}

#[no_mangle]
pub unsafe extern "C" fn setTimeCourseSelectionList(handle: RRHandle, list: *const c_char) -> bool {
    set_selections(handle, false, || Ok(split_selections(&st(list)?)))
}

#[no_mangle]
pub unsafe extern "C" fn setTimeCourseSelectionListEx(handle: RRHandle, length: c_int, list: *mut *const c_char) -> bool {
    set_selections(handle, false, || strs(list, length))
}

#[no_mangle]
pub unsafe extern "C" fn getTimeCourseSelectionList(handle: RRHandle) -> *mut RRStringArray {
    with(handle, null_mut(), |r| Ok(list_or_null(r.rr.selections())))
}

#[no_mangle]
pub unsafe extern "C" fn setSteadyStateSelectionList(handle: RRHandle, list: *const c_char) -> bool {
    set_selections(handle, true, || Ok(split_selections(&st(list)?)))
}

#[no_mangle]
pub unsafe extern "C" fn setSteadyStateSelectionListEx(handle: RRHandle, length: c_int, list: *mut *const c_char) -> bool {
    set_selections(handle, true, || strs(list, length))
}

#[no_mangle]
pub unsafe extern "C" fn getSteadyStateSelectionList(handle: RRHandle) -> *mut RRStringArray {
    with(handle, null_mut(), |r| Ok(string_array(r.rr.steady_state_selections())))
}

// ---------------------------------------------------------------- simulation

#[no_mangle]
pub unsafe extern "C" fn simulate(handle: RRHandle) -> *mut RRCData {
    wm(handle, null_mut(), |rr| Ok(cdata_of(rr.simulate()?, None)))
}

#[no_mangle]
pub unsafe extern "C" fn simulateNoReturn(handle: RRHandle) -> bool {
    wm(handle, false, |rr| { rr.simulate()?; Ok(true) })
}

/// The last result (no rows before the first simulation).
#[no_mangle]
pub unsafe extern "C" fn getSimulationResult(handle: RRHandle) -> *mut RRCData {
    with(handle, null_mut(), |r| Ok(match r.rr.result() {
        Some(m) => cdata_of(m, None),
        None => cdata::<String>(0, &[], None, r.rr.selections()),
    }))
}

/// Upstream returns its C++ ls::DoubleMatrix; here an RRDoubleMatrix* owned by the instance (do not
/// free; valid until the next call or freeRRInstance). NULL before the first simulation.
#[no_mangle]
pub unsafe extern "C" fn getSimulationResultAsDoubleMatrix(handle: RRHandle) -> RRHandle {
    with(handle, null_mut(), |r| {
        crate::utilities::freeMatrix(r.last);
        r.last = r.rr.result().map_or(null_mut(), |m| matrix(m.nr, m.nc, &m.data));
        Ok(r.last as RRHandle)
    })
}

unsafe fn grid(handle: RRHandle, start: f64, end: f64, points: c_int) -> bool {
    setTimeStart(handle, start) && setTimeEnd(handle, end) && setNumPoints(handle, points)
}

#[no_mangle]
pub unsafe extern "C" fn simulateEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> *mut RRCData {
    if grid(handle, timeStart, timeEnd, numberOfPoints) { simulate(handle) } else { null_mut() }
}

#[no_mangle]
pub unsafe extern "C" fn simulateExNoReturn(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> bool {
    grid(handle, timeStart, timeEnd, numberOfPoints) && simulateNoReturn(handle)
}

#[no_mangle]
pub unsafe extern "C" fn simulateTimes(handle: RRHandle, times: *const f64, size: c_int) -> *mut RRCData {
    if setTimes(handle, times, size) { simulate(handle) } else { null_mut() }
}

#[no_mangle]
pub unsafe extern "C" fn oneStep(handle: RRHandle, currentTime: f64, stepSize: f64, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.one_step(currentTime, stepSize)?))
}

// ---------------------------------------------------------------- steady state

#[no_mangle]
pub unsafe extern "C" fn steadyState(handle: RRHandle, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.steady_state()?))
}

#[no_mangle]
pub unsafe extern "C" fn computeSteadyStateValues(handle: RRHandle) -> *mut RRVector {
    wm(handle, null_mut(), |rr| Ok(vector(&rr.steady_state_values()?)))
}

// ---------------------------------------------------------------- get / set

#[no_mangle]
pub unsafe extern "C" fn getValue(handle: RRHandle, symbolId: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.get_value(&st(symbolId)?)?))
}

#[no_mangle]
pub unsafe extern "C" fn setValue(handle: RRHandle, symbolId: *const c_char, value: f64) -> bool {
    wm(handle, false, |rr| { rr.set_value(&st(symbolId)?, value)?; Ok(true) })
}

macro_rules! vectors {
    ($($name:ident => |$rr:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle) -> *mut RRVector { wm(handle, null_mut(), |$rr| Ok(vector(&$body))) }
    )*};
}

/// Floating species values `first..first+n` in the current order (upstream: the independent species
/// are the first getNumberOfIndependentSpecies ones, the dependent ones follow).
fn floating_part(rr: &mut RoadRunner, amounts: bool, dependent: bool) -> R<Vec<f64>> {
    let all = if amounts { rr.floating_amounts()? } else { rr.floating_concentrations()? };
    let (ni, nd) = (rr.independent_ids().len(), rr.dependent_ids().len());
    Ok(if dependent { all[ni..ni + nd].to_vec() } else { all[..ni].to_vec() })
}

vectors! {
    getFloatingSpeciesConcentrations => |rr| rr.floating_concentrations()?;
    getFloatingSpeciesAmounts => |rr| rr.floating_amounts()?;
    getRateRuleValues => |rr| rr.rate_rule_values()?;
    getIndependentFloatingSpeciesConcentrations => |rr| floating_part(rr, false, false)?;
    getDependentFloatingSpeciesConcentrations => |rr| floating_part(rr, false, true)?;
    getIndependentFloatingSpeciesAmounts => |rr| floating_part(rr, true, false)?;
    getDependentFloatingSpeciesAmounts => |rr| floating_part(rr, true, true)?;
    getBoundarySpeciesConcentrations => |rr| rr.boundary_concentrations()?;
    getBoundarySpeciesAmounts => |rr| rr.boundary_amounts()?;
    getGlobalParameterValues => |rr| rr.global_parameter_values()?;
    getFloatingSpeciesInitialConcentrations => |rr| rr.floating_init_concentrations()?;
    getReactionRates => |rr| rr.reaction_rates()?;
}

/// Rates of change of getRatesOfChangeIds (rate-rule variables, then the floating species amounts;
/// upstream returns only the independent floating species', misaligned with its ids). NULL if none.
#[no_mangle]
pub unsafe extern "C" fn getRatesOfChange(handle: RRHandle) -> *mut RRVector {
    wm(handle, null_mut(), |rr| { let v = rr.rates_of_change()?; Ok(if v.is_empty() { null_mut() } else { vector(&v) }) })
}

/// The value of selection `fmt(ids[index])`.
fn by_index(rr: &mut RoadRunner, ids: Vec<String>, index: c_int, f: &str, what: &str) -> R<f64> {
    let i = idx(index, ids.len(), what)?;
    rr.get_value(&f.replace("{}", &ids[i]))
}

fn set_by_index(rr: &mut RoadRunner, ids: Vec<String>, index: c_int, f: &str, what: &str, v: f64) -> R<bool> {
    let i = idx(index, ids.len(), what)?;
    rr.set_value(&f.replace("{}", &ids[i]), v)?;
    Ok(true)
}

#[no_mangle]
pub unsafe extern "C" fn setBoundarySpeciesByIndex(handle: RRHandle, index: c_uint, value: f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.boundary_ids(); set_by_index(rr, ids, index.min(i32::MAX as u32) as c_int, "[{}]", "setBoundarySpeciesByIndex", value) })
}

#[no_mangle]
pub unsafe extern "C" fn getBoundarySpeciesByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.boundary_ids(); out(value, by_index(rr, ids, index, "[{}]", "getBoundarySpeciesByIndex")?) })
}

/// Initial concentration of floating species `index`, then reset() (parameters keep their values).
#[no_mangle]
pub unsafe extern "C" fn setFloatingSpeciesInitialConcentrationByIndex(handle: RRHandle, index: c_int, value: f64) -> bool {
    wm(handle, false, |rr| {
        let ids = rr.floating_ids();
        rr.set_init(&ids[idx(index, ids.len(), "setFloatingSpeciesInitialConcentrationByIndex")?], true, value)?;
        rr.reset()?;
        Ok(true)
    })
}

/// Out of range is an error (upstream returns true without writing *value).
#[no_mangle]
pub unsafe extern "C" fn getFloatingSpeciesInitialConcentrationByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.floating_ids(); out(value, by_index(rr, ids, index, "init([{}])", "getFloatingSpeciesInitialConcentrationByIndex")?) })
}

#[no_mangle]
pub unsafe extern "C" fn setFloatingSpeciesByIndex(handle: RRHandle, index: c_int, value: f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.floating_ids(); set_by_index(rr, ids, index, "[{}]", "setFloatingSpeciesByIndex", value) })
}

#[no_mangle]
pub unsafe extern "C" fn getFloatingSpeciesByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.floating_ids(); out(value, by_index(rr, ids, index, "[{}]", "getFloatingSpeciesByIndex")?) })
}

#[no_mangle]
pub unsafe extern "C" fn setGlobalParameterByIndex(handle: RRHandle, index: c_int, value: f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.global_parameter_ids(); set_by_index(rr, ids, index, "{}", "setGlobalParameterByIndex", value) })
}

#[no_mangle]
pub unsafe extern "C" fn getGlobalParameterByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.global_parameter_ids(); out(value, by_index(rr, ids, index, "{}", "getGlobalParameterByIndex")?) })
}

#[no_mangle]
pub unsafe extern "C" fn getCompartmentByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.compartment_ids(); out(value, by_index(rr, ids, index, "{}", "getCompartmentByIndex")?) })
}

#[no_mangle]
pub unsafe extern "C" fn setCompartmentByIndex(handle: RRHandle, index: c_int, value: f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.compartment_ids(); set_by_index(rr, ids, index, "{}", "setCompartmentByIndex", value) })
}

#[no_mangle]
pub unsafe extern "C" fn getReactionRate(handle: RRHandle, index: c_int, rate: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.reaction_ids(); out(rate, by_index(rr, ids, index, "{}", "getReactionRate")?) })
}

/// Amount rate of change of floating species `index` (as upstream).
#[no_mangle]
pub unsafe extern "C" fn getRateOfChange(handle: RRHandle, index: c_int, value: *mut f64) -> bool {
    wm(handle, false, |rr| { let ids = rr.floating_ids(); out(value, by_index(rr, ids, index, "{}'", "getRateOfChange")?) })
}

/// Set the first len(vec) values of `ids` (selection format `f`).
fn set_all(rr: &mut RoadRunner, ids: Vec<String>, v: &[f64], f: &str) -> R<()> {
    if v.len() > ids.len() { return Err(format!("{} values for {} species", v.len(), ids.len())); }
    ids.iter().zip(v).try_for_each(|(s, &x)| rr.set_value(&f.replace("{}", s), x))
}

#[no_mangle]
pub unsafe extern "C" fn setFloatingSpeciesConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool {
    wm(handle, false, |rr| { let ids = rr.floating_ids(); set_all(rr, ids, &vec_arg(vec)?, "[{}]")?; Ok(true) })
}

#[no_mangle]
pub unsafe extern "C" fn setBoundarySpeciesConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool {
    wm(handle, false, |rr| { let ids = rr.boundary_ids(); set_all(rr, ids, &vec_arg(vec)?, "[{}]")?; Ok(true) })
}

/// changeInitialConditions: initial and current concentrations of the first len(vec) floating
/// species; nothing is reset.
#[no_mangle]
pub unsafe extern "C" fn setFloatingSpeciesInitialConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool {
    wm(handle, false, |rr| {
        let (ids, v) = (rr.floating_ids(), vec_arg(vec)?);
        if v.len() > ids.len() { return Err(format!("{} values for {} species", v.len(), ids.len())); }
        for (s, &x) in ids.iter().zip(&v) { rr.set_init(s, true, x)?; }
        set_all(rr, ids, &v, "[{}]")?;
        Ok(true)
    })
}

/// Sets the floating species concentrations (they stay), returns the reaction rates.
#[no_mangle]
pub unsafe extern "C" fn getReactionRatesEx(handle: RRHandle, vec: *mut RRVector) -> *mut RRVector {
    wm(handle, null_mut(), |rr| { let ids = rr.floating_ids(); set_all(rr, ids, &vec_arg(vec)?, "[{}]")?; Ok(vector(&rr.reaction_rates()?)) })
}

/// Sets the floating species concentrations (they stay), returns getRatesOfChange.
#[no_mangle]
pub unsafe extern "C" fn getRatesOfChangeEx(handle: RRHandle, vec: *mut RRVector) -> *mut RRVector {
    wm(handle, null_mut(), |rr| { let ids = rr.floating_ids(); set_all(rr, ids, &vec_arg(vec)?, "[{}]")?; Ok(vector(&rr.rates_of_change()?)) })
}

/// The engine evaluates on demand: only checks that a model is loaded.
#[no_mangle]
pub unsafe extern "C" fn evalModel(handle: RRHandle) -> bool { wm(handle, false, |_| Ok(true)) }

// ---------------------------------------------------------------- matrices

macro_rules! matrices {
    ($($name:ident => |$rr:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle) -> *mut RRDoubleMatrix { wm(handle, null_mut(), |$rr| Ok(mat($body))) }
    )*};
}

matrices! {
    getFullJacobian => |rr| rr.full_jacobian()?;
    getReducedJacobian => |rr| rr.reduced_jacobian()?;
    getStoichiometryMatrix => |rr| rr.stoichiometry_matrix()?;
    getLinkMatrix => |rr| rr.link_matrix()?;
    getNrMatrix => |rr| rr.reduced_stoichiometry_matrix()?;
    getUnscaledElasticityMatrix => |rr| rr.elasticity_matrix(false)?;
    getScaledElasticityMatrix => |rr| rr.elasticity_matrix(true)?;
    getUnscaledConcentrationControlCoefficientMatrix => |rr| rr.cc_matrix(false)?;
    getScaledConcentrationControlCoefficientMatrix => |rr| rr.cc_matrix(true)?;
    getUnscaledFluxControlCoefficientMatrix => |rr| rr.fcc_matrix(false)?;
    getScaledFluxControlCoefficientMatrix => |rr| rr.fcc_matrix(true)?;
}

/// Gamma (moiety totals x floating species); not NULL when empty (as upstream).
#[no_mangle]
pub unsafe extern "C" fn getConservationMatrix(handle: RRHandle) -> *mut RRDoubleMatrix {
    wm(handle, null_mut(), |rr| { let m = rr.conservation_matrix()?; Ok(matrix(m.nr, m.nc, &m.data)) })
}

/// Eigenvalues of the full Jacobian: n x 2 (real, imaginary); NULL when there are none.
#[no_mangle]
pub unsafe extern "C" fn getEigenvalues(handle: RRHandle) -> *mut RRDoubleMatrix {
    wm(handle, null_mut(), |rr| {
        let e = rr.full_eigenvalues()?;
        Ok(if e.is_empty() { null_mut() } else { matrix(e.len(), 2, &e.iter().flat_map(|z| [z.0, z.1]).collect::<Vec<_>>()) })
    })
}

// ---------------------------------------------------------------- counts and ids

macro_rules! counts {
    ($($name:ident => |$rr:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle) -> c_int { wm(handle, -1, |$rr| Ok($body.len() as c_int)) }
    )*};
}

counts! {
    getNumberOfReactions => |rr| rr.reaction_ids();
    getNumberOfCompartments => |rr| rr.compartment_ids();
    getNumberOfBoundarySpecies => |rr| rr.boundary_ids();
    getNumberOfFloatingSpecies => |rr| rr.floating_ids();
    getNumberOfGlobalParameters => |rr| rr.global_parameter_ids();
    getNumberOfDependentSpecies => |rr| rr.dependent_ids();
    getNumberOfIndependentSpecies => |rr| rr.independent_ids();
}

macro_rules! id_lists {
    ($($name:ident $mk:ident => |$rr:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle) -> *mut RRStringArray { wm(handle, null_mut(), |$rr| Ok($mk(&$body))) }
    )*};
}

fn always(v: &[String]) -> *mut RRStringArray { string_array(v) }

id_lists! {
    getReactionIds list_or_null => |rr| rr.reaction_ids();
    getBoundarySpeciesIds list_or_null => |rr| rr.boundary_ids();
    getBoundarySpeciesConcentrationIds list_or_null => |rr| conc(rr.boundary_ids());
    getFloatingSpeciesIds list_or_null => |rr| rr.floating_ids();
    getDependentFloatingSpeciesIds list_or_null => |rr| rr.dependent_ids();
    getIndependentFloatingSpeciesIds list_or_null => |rr| rr.independent_ids();
    getFloatingSpeciesConcentrationIds list_or_null => |rr| conc(rr.floating_ids());
    getGlobalParameterIds list_or_null => |rr| rr.global_parameter_ids();
    getRatesOfChangeIds list_or_null => |rr| rr.rates_of_change_ids();
    getCompartmentIds always => |rr| rr.compartment_ids();
    getEigenvalueIds always => |rr| rr.eigenvalue_ids();
    getFloatingSpeciesInitialConditionIds always => |rr| rr.floating_init_ids();
    getFloatingSpeciesInitialConcentrationIds always => |rr| rr.floating_init_ids();
}

macro_rules! rr_lists {
    ($($name:ident => |$rr:ident| $body:expr;)*) => {$(
        #[no_mangle]
        pub unsafe extern "C" fn $name(handle: RRHandle) -> *mut RRList { wm(handle, null_mut(), |$rr| Ok(rrlist(&$body))) }
    )*};
}

/// Upstream's group list ({label, {items}} per group; nested coefficient groups {label, {{name, {items}}, ...}}).
fn symbols(rr: &RoadRunner, steady: bool) -> Vec<L> {
    let (fs, bs) = (rr.floating_ids(), rr.boundary_ids());
    let nested = |label: &str, v| L::L(vec![L::S(label.into()), L::L(groups(v))]);
    let mut v = vec![];
    if !steady { v.push(group("Time", vec!["time".into()])); }
    v.push(group("Floating Species", conc(fs.clone())));
    v.push(group("Boundary Species", conc(bs.clone())));
    v.push(group("Floating Species (amount)", fs));
    v.push(group("Boundary Species (amount)", bs));
    v.push(group("Global Parameters", rr.global_parameter_ids()));
    if steady {
        v.push(group("Volumes", rr.compartment_ids()));
        v.push(group("Fluxes", rr.reaction_ids()));
        v.push(nested("Flux Control Coefficients", rr.cc_ids(true, true)));
        v.push(nested("Concentration Control Coefficients", rr.cc_ids(false, true)));
        v.push(nested("Unscaled Concentration Control Coefficients", rr.cc_ids(false, false)));
    } else {
        v.push(group("Fluxes", rr.reaction_ids()));
        v.push(group("Rates of Change", rr.rates_of_change_ids()));
        v.push(group("Volumes", rr.compartment_ids()));
    }
    v.push(nested("Elasticity Coefficients", rr.ec_ids(true)));
    v.push(nested("Unscaled Elasticity Coefficients", rr.ec_ids(false)));
    v.push(group("Eigenvalues", rr.eigenvalue_ids()));
    v
}

rr_lists! {
    getAvailableTimeCourseSymbols => |rr| symbols(rr, false);
    getAvailableSteadyStateSymbols => |rr| symbols(rr, true);
    getElasticityCoefficientIds => |rr| groups(rr.ec_ids(true));
    getUnscaledElasticityCoefficientIds => |rr| groups(rr.ec_ids(false));
    getFluxControlCoefficientIds => |rr| groups(rr.cc_ids(true, true));
    getUnscaledFluxControlCoefficientIds => |rr| groups(rr.cc_ids(true, false));
    getConcentrationControlCoefficientIds => |rr| groups(rr.cc_ids(false, true));
    getUnscaledConcentrationControlCoefficientIds => |rr| groups(rr.cc_ids(false, false));
}

// ---------------------------------------------------------------- MCA

#[no_mangle]
pub unsafe extern "C" fn getScaledFloatingSpeciesElasticity(handle: RRHandle, reactionId: *const c_char, speciesId: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| {
        let s = st(speciesId)?;
        if !rr.floating_ids().contains(&s) { return Err(format!("invalid species name: {s}")); }
        out(value, rr.get_ee(&st(reactionId)?, &s, false)?)
    })
}

#[no_mangle]
pub unsafe extern "C" fn getuCC(handle: RRHandle, variable: *const c_char, parameter: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.get_ucc(&st(variable)?, &st(parameter)?)?))
}

#[no_mangle]
pub unsafe extern "C" fn getCC(handle: RRHandle, variable: *const c_char, parameter: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.get_cc(&st(variable)?, &st(parameter)?)?))
}

/// At steady state (upstream's default).
#[no_mangle]
pub unsafe extern "C" fn getEE(handle: RRHandle, name: *const c_char, species: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.get_ee(&st(name)?, &st(species)?, true)?))
}

#[no_mangle]
pub unsafe extern "C" fn getuEE(handle: RRHandle, name: *const c_char, species: *const c_char, value: *mut f64) -> bool {
    wm(handle, false, |rr| out(value, rr.get_uee(&st(name)?, &st(species)?, true)?))
}

#[no_mangle]
pub unsafe extern "C" fn getDiffStepSize(handle: RRHandle, value: *mut f64) -> bool { with(handle, false, |r| out(value, r.rr.diff_step)) }

#[no_mangle]
pub unsafe extern "C" fn setDiffStepSize(handle: RRHandle, value: f64) -> bool {
    with(handle, false, |r| { r.rr.diff_step = value; Ok(true) })
}

// ---------------------------------------------------------------- stochastic

#[no_mangle]
pub unsafe extern "C" fn getSeed(handle: RRHandle, seed: *mut i64) -> bool { with(handle, false, |r| out(seed, r.rr.seed() as i64)) }

#[no_mangle]
pub unsafe extern "C" fn setSeed(handle: RRHandle, seed: i64) -> bool {
    with(handle, false, |r| { r.rr.set_seed(seed as u64)?; Ok(true) })
}

/// Make gillespie current with the given variable_step_size (as upstream, it stays current).
fn use_gillespie(r: &mut Inst, variable: bool) -> R<()> {
    r.set_integrator("gillespie")?;
    r.rr.integrator_mut().set("variable_step_size", Value::Bool(variable))
}

/// Variable steps: one row per reaction firing (upstream sets variable_step_size = false here,
/// contradicting its documentation, so its gillespie() equals gillespieOnGrid()).
#[no_mangle]
pub unsafe extern "C" fn gillespie(handle: RRHandle) -> *mut RRCData {
    with(handle, null_mut(), |r| { use_gillespie(r, true)?; Ok(cdata_of(r.rr.simulate()?, None)) })
}

#[no_mangle]
pub unsafe extern "C" fn gillespieEx(handle: RRHandle, timeStart: f64, timeEnd: f64) -> *mut RRCData {
    if setTimeStart(handle, timeStart) && setTimeEnd(handle, timeEnd) { gillespie(handle) } else { null_mut() }
}

#[no_mangle]
pub unsafe extern "C" fn gillespieOnGrid(handle: RRHandle) -> *mut RRCData {
    with(handle, null_mut(), |r| { use_gillespie(r, false)?; Ok(cdata_of(r.rr.simulate()?, None)) })
}

#[no_mangle]
pub unsafe extern "C" fn gillespieOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> *mut RRCData {
    if grid(handle, timeStart, timeEnd, numberOfPoints) { gillespieOnGrid(handle) } else { null_mut() }
}

/// n runs from the reset state on the grid: (mean, sample sd); zeros for n <= 0 (as upstream).
fn gillespie_stats(r: &mut Inst, n: c_int) -> R<(Matrix, Vec<f64>)> {
    use_gillespie(r, false)?;
    let rr = &mut r.rr;
    if !rr.is_model_loaded() { return Err(NO_MODEL.into()); }
    rr.reset()?;
    let (nr, nc) = (rr.points, rr.selections().len());
    let per = nr * nc;
    let (mut mean, mut sd) = (vec![0.0; per], vec![0.0; per]);
    if n > 0 {
        let n = n as usize;
        let y = rr.gillespie_runs(n, rr.start, rr.end, nr, None)?;
        for run in y.chunks(per) { for (m, v) in mean.iter_mut().zip(run) { *m += v / n as f64; } }
        if n > 1 {
            for run in y.chunks(per) { for k in 0..per { sd[k] += (run[k] - mean[k]).powi(2) / (n - 1) as f64; } }
            sd.iter_mut().for_each(|v| *v = v.sqrt());
        }
    }
    Ok((Matrix { nr, nc, data: mean, rows: vec![], cols: rr.selections().to_vec() }, sd))
}

#[no_mangle]
pub unsafe extern "C" fn gillespieMeanOnGrid(handle: RRHandle, numberOfSimulations: c_int) -> *mut RRCData {
    with(handle, null_mut(), |r| Ok(cdata_of(&gillespie_stats(r, numberOfSimulations)?.0, None)))
}

#[no_mangle]
pub unsafe extern "C" fn gillespieMeanOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int, numberOfSimulations: c_int) -> *mut RRCData {
    if grid(handle, timeStart, timeEnd, numberOfPoints) { gillespieMeanOnGrid(handle, numberOfSimulations) } else { null_mut() }
}

/// Data: mean, Weights: sample standard deviation.
#[no_mangle]
pub unsafe extern "C" fn gillespieMeanSDOnGrid(handle: RRHandle, numberOfSimulations: c_int) -> *mut RRCData {
    with(handle, null_mut(), |r| { let (m, sd) = gillespie_stats(r, numberOfSimulations)?; Ok(cdata_of(&m, Some(&sd))) })
}

#[no_mangle]
pub unsafe extern "C" fn gillespieMeanSDOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfSteps: c_int, numberOfSimulations: c_int) -> *mut RRCData {
    if grid(handle, timeStart, timeEnd, numberOfSteps) { gillespieMeanSDOnGrid(handle, numberOfSimulations) } else { null_mut() }
}

// ---------------------------------------------------------------- reset

#[no_mangle]
pub unsafe extern "C" fn reset(handle: RRHandle) -> bool { wm(handle, false, |rr| { rr.reset()?; Ok(true) }) }

#[no_mangle]
pub unsafe extern "C" fn resetAll(handle: RRHandle) -> bool { wm(handle, false, |rr| { rr.reset_all()?; Ok(true) }) }

#[no_mangle]
pub unsafe extern "C" fn resetToOrigin(handle: RRHandle) -> bool { wm(handle, false, |rr| { rr.reset_to_origin()?; Ok(true) }) }

#[no_mangle]
pub unsafe extern "C" fn resetParameter(handle: RRHandle) -> bool { wm(handle, false, |rr| { rr.reset_parameter()?; Ok(true) }) }

// ---------------------------------------------------------------- ExecutableModel access (RRPlugins section)

#[no_mangle]
pub unsafe extern "C" fn _getNumIndFloatingSpecies(handle: RRHandle) -> c_int {
    wm(handle, -1, |rr| Ok((if rr.conservation_laws() { rr.independent_ids() } else { rr.floating_ids() }).len() as c_int))
}

#[no_mangle]
pub unsafe extern "C" fn _getNumRateRules(handle: RRHandle) -> c_int { wm(handle, -1, |rr| Ok(rr.rate_rule_ids().len() as c_int)) }

#[no_mangle]
pub unsafe extern "C" fn _getTime(handle: RRHandle) -> f64 { wm(handle, -1.0, |rr| Ok(rr.time())) }

/// Size of the state vector (rate-rule variables and the independent floating species).
#[no_mangle]
pub unsafe extern "C" fn _getStateVector(handle: RRHandle) -> c_int { wm(handle, -1, |rr| Ok(rr.rates_of_change_ids().len() as c_int)) }

/// Rates of the state vector at the current state and time `time` into `value` (_getStateVector() entries).
#[no_mangle]
pub unsafe extern "C" fn _getStateVectorRate(handle: RRHandle, time: f64, value: *mut f64) {
    wm(handle, (), |rr| {
        if value.is_null() { return Err("NULL output argument".into()); }
        let t0 = rr.time();
        rr.set_value("time", time)?;
        let v = rr.rates_of_change();
        rr.set_value("time", t0)?;
        let v = v?;
        value.copy_from_nonoverlapping(v.as_ptr(), v.len());
        Ok(())
    })
}

// ---------------------------------------------------------------- writeRRData (rrc_utilities.h)

/// The last simulation result as CSV (header: the selections).
pub(crate) unsafe fn write_rr_data(handle: RRHandle, path: *const c_char) -> c_int {
    with(handle, 0, |r| {
        let m = r.rr.result().ok_or("no simulation result")?;
        let mut s = m.cols.join(",") + "\n";
        for row in m.data.chunks(m.nc.max(1)) { s += &(row.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(",") + "\n"); }
        let p = st(path)?;
        std::fs::write(&p, s).map_err(|e| format!("{p}: {e}"))?;
        Ok(1)
    })
}
