//! rrc_api.h. Instance-independent functions are implemented here; everything that needs a
//! RoadRunner instance is a stub (bottom of file) that records "not implemented yet" and returns
//! the upstream failure value (false / NULL / -1). Replace a stub by deleting its line and writing the fn.
use crate::types::*;
use crate::{has_error, last_error_ptr, set_error};
use libc::{c_char, c_int, c_uint};
use std::collections::BTreeMap;
use std::ptr::null_mut;
use std::sync::Mutex;

pub(crate) fn not_implemented<T>(name: &str, fail: T) -> T {
    set_error(format!("{name}: not implemented yet"));
    fail
}

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
pub extern "C" fn getInfo(_handle: RRHandle) -> *mut c_char {
    text(&format!("roadrunner_nitro C API {}\nTemp folder: {}\n", version_ex(), folder("temp")))
}

#[no_mangle]
pub extern "C" fn getlibSBMLVersion(_handle: RRHandle) -> *mut c_char {
    text("none")
}

// ---------------------------------------------------------------- folders & compiler (process-wide, stored only)
// ponytail: one process-wide store ignoring the handle; move into the instance when the core lands.

static FOLDERS: Mutex<BTreeMap<&'static str, String>> = Mutex::new(BTreeMap::new());

fn folder(key: &'static str) -> String {
    FOLDERS.lock().unwrap().get(key).cloned().unwrap_or_else(|| match key {
        "temp" => std::env::temp_dir().to_string_lossy().into_owned(),
        _ => String::new(),
    })
}

unsafe fn set_folder(key: &'static str, value: *const c_char) -> bool {
    let Some(v) = cstr(value) else {
        set_error(format!("{key}: NULL argument"));
        return false;
    };
    FOLDERS.lock().unwrap().insert(key, v);
    true
}

#[no_mangle]
pub unsafe extern "C" fn setInstallFolder(folder: *const c_char) -> bool {
    set_folder("install", folder)
}

#[no_mangle]
pub unsafe extern "C" fn setTempFolder(_handle: RRHandle, folder: *const c_char) -> bool {
    set_folder("temp", folder)
}

#[no_mangle]
pub extern "C" fn getTempFolder(_handle: RRHandle) -> *mut c_char {
    text(&folder("temp"))
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

#[no_mangle]
pub unsafe extern "C" fn setCompiler(_handle: RRHandle, fNameWithPath: *const c_char) -> bool {
    set_folder("compiler", fNameWithPath)
}

#[no_mangle]
pub extern "C" fn getCompiler(_handle: RRHandle) -> *mut c_char {
    text(&folder("compiler"))
}

#[no_mangle]
pub unsafe extern "C" fn setCompilerLocation(_handle: RRHandle, folder: *const c_char) -> bool {
    set_folder("compilerLocation", folder)
}

#[no_mangle]
pub extern "C" fn getCompilerLocation(_handle: RRHandle) -> *mut c_char {
    text(&folder("compilerLocation"))
}

#[no_mangle]
pub unsafe extern "C" fn setSupportCodeFolder(_handle: RRHandle, folder: *const c_char) -> bool {
    set_folder("supportCodeFolder", folder)
}

#[no_mangle]
pub extern "C" fn getSupportCodeFolder(_handle: RRHandle) -> *mut c_char {
    text(&folder("supportCodeFolder"))
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

// ---------------------------------------------------------------- stubs (need a RoadRunner instance)

macro_rules! stubs {
    ($($name:ident($($arg:ident: $t:ty),*) -> $r:ty = $fail:expr;)*) => {$(
        #[no_mangle]
        pub extern "C" fn $name($($arg: $t),*) -> $r { $(let _ = $arg;)* not_implemented(stringify!($name), $fail) }
    )*};
}

stubs! {
    // Initialization, version info, folders, compiler, errors
    createRRInstance() -> RRHandle = null_mut();
    createRRInstanceEx(tempFolder: *const c_char, compiler: *const c_char) -> RRHandle = null_mut();
    freeRRInstance(handle: RRHandle) -> bool = false;
    // Flags/Options
    setComputeAndAssignConservationLaws(handle: RRHandle, On_Or_Off: bool) -> bool = false;
    getComputeAndAssignConservationLaws(handle: RRHandle, value: *mut c_int) -> c_int = -1;
    // Read and Write models
    loadSBML(handle: RRHandle, sbml: *const c_char) -> bool = false;
    loadSBMLEx(handle: RRHandle, sbml: *const c_char, forceRecompile: bool) -> bool = false;
    saveState(handle: RRHandle, filename: *const c_char) -> bool = false;
    loadState(handle: RRHandle, filename: *const c_char) -> bool = false;
    loadSBMLFromFile(handle: RRHandle, fileName: *const c_char) -> bool = false;
    loadSBMLFromFileE(handle: RRHandle, fileName: *const c_char, forceRecompile: bool) -> bool = false;
    clearModel(handle: RRHandle) -> bool = false;
    isModelLoaded(handle: RRHandle) -> bool = false;
    loadSimulationSettings(handle: RRHandle, fileName: *const c_char) -> bool = false;
    getCurrentSBML(handle: RRHandle) -> *mut c_char = null_mut();
    getSBML(handle: RRHandle) -> *mut c_char = null_mut();
    // Model editing methods
    addSpeciesAmount(handle: RRHandle, sid: *const c_char, compartment: *const c_char, initialAmount: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool) -> bool = false;
    addSpeciesConcentration(handle: RRHandle, sid: *const c_char, compartment: *const c_char, initialConcentration: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool) -> bool = false;
    addSpeciesConcentrationNoRegen(handle: RRHandle, sid: *const c_char, compartment: *const c_char, initialConcentration: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool) -> bool = false;
    addSpeciesAmountNoRegen(handle: RRHandle, sid: *const c_char, compartment: *const c_char, initialAmount: f64, hasOnlySubstanceUnits: bool, boundaryCondition: bool) -> bool = false;
    removeSpecies(handle: RRHandle, sid: *const c_char) -> bool = false;
    removeSpeciesNoRegen(handle: RRHandle, sid: *const c_char) -> bool = false;
    setBoundary(handle: RRHandle, sid: *const c_char, boundaryCondition: bool) -> bool = false;
    setBoundaryNoRegen(handle: RRHandle, sid: *const c_char, boundaryCondition: bool) -> bool = false;
    setHasOnlySubstanceUnits(handle: RRHandle, sid: *const c_char, hasOnlySubstanceUnits: bool) -> bool = false;
    setHasOnlySubstanceUnitsNoRegen(handle: RRHandle, sid: *const c_char, hasOnlySubstanceUnits: bool) -> bool = false;
    setInitAmount(handle: RRHandle, sid: *const c_char, initAmount: f64) -> bool = false;
    setInitConcentration(handle: RRHandle, sid: *const c_char, initConcentration: f64) -> bool = false;
    setConstant(handle: RRHandle, sid: *const c_char, constant: bool) -> bool = false;
    setConstantNoRegen(handle: RRHandle, sid: *const c_char, constant: bool) -> bool = false;
    addReactionFromSBML(handle: RRHandle, sbmlRep: *const c_char) -> bool = false;
    addReactionFromSBMLNoRegen(handle: RRHandle, sbmlRep: *const c_char) -> bool = false;
    addReaction(handle: RRHandle, rid: *const c_char, reactants: *mut *const c_char, numReactants: c_int, products: *mut *const c_char, numProducts: c_int, kineticLaw: *const c_char) -> bool = false;
    addReactionNoRegen(handle: RRHandle, rid: *const c_char, reactants: *mut *const c_char, numReactants: c_int, products: *mut *const c_char, numProducts: c_int, kineticLaw: *const c_char) -> bool = false;
    removeReaction(handle: RRHandle, rid: *const c_char) -> bool = false;
    removeReactionNoRegen(handle: RRHandle, rid: *const c_char) -> bool = false;
    setReversible(handle: RRHandle, rid: *const c_char, reversible: bool) -> bool = false;
    setReversibleNoRegen(handle: RRHandle, rid: *const c_char, reversible: bool) -> bool = false;
    setKineticLaw(handle: RRHandle, rid: *const c_char, kineticLaw: *const c_char) -> bool = false;
    setKineticLawNoRegen(handle: RRHandle, rid: *const c_char, kineticLaw: *const c_char) -> bool = false;
    getKineticLaw(handle: RRHandle, rid: *const c_char) -> *mut c_char = null_mut();
    addParameter(handle: RRHandle, pid: *const c_char, value: f64) -> bool = false;
    addParameterNoRegen(handle: RRHandle, pid: *const c_char, value: f64) -> bool = false;
    removeParameter(handle: RRHandle, pid: *const c_char) -> bool = false;
    removeParameterNoRegen(handle: RRHandle, pid: *const c_char) -> bool = false;
    addCompartment(handle: RRHandle, cid: *const c_char, initVolume: f64) -> bool = false;
    addCompartmentNoRegen(handle: RRHandle, cid: *const c_char, initVolume: f64) -> bool = false;
    removeCompartment(handle: RRHandle, cid: *const c_char) -> bool = false;
    removeCompartmentNoRegen(handle: RRHandle, cid: *const c_char) -> bool = false;
    addInitialAssignment(handle: RRHandle, vid: *const c_char, formula: *const c_char, forceRegenerate: bool) -> bool = false;
    addAssignmentRule(handle: RRHandle, vid: *const c_char, formula: *const c_char) -> bool = false;
    addAssignmentRuleNoRegen(handle: RRHandle, vid: *const c_char, formula: *const c_char) -> bool = false;
    addRateRule(handle: RRHandle, vid: *const c_char, formula: *const c_char) -> bool = false;
    addRateRuleNoRegen(handle: RRHandle, vid: *const c_char, formula: *const c_char) -> bool = false;
    removeRules(handle: RRHandle, vid: *const c_char) -> bool = false;
    removeRulesNoRegen(handle: RRHandle, vid: *const c_char) -> bool = false;
    addEvent(handle: RRHandle, eid: *const c_char, useValuesFromTriggerTime: bool, trigger: *const c_char) -> bool = false;
    addEventNoRegen(handle: RRHandle, eid: *const c_char, useValuesFromTriggerTime: bool, trigger: *const c_char) -> bool = false;
    addTrigger(handle: RRHandle, eid: *const c_char, trigger: *const c_char) -> bool = false;
    addTriggerNoRegen(handle: RRHandle, eid: *const c_char, trigger: *const c_char) -> bool = false;
    setPersistent(handle: RRHandle, eid: *const c_char, persistent: bool) -> bool = false;
    setPersistentNoRegen(handle: RRHandle, eid: *const c_char, persistent: bool) -> bool = false;
    setTriggerInitialValue(handle: RRHandle, eid: *const c_char, initValue: bool) -> bool = false;
    setTriggerInitialValueNoRegen(handle: RRHandle, eid: *const c_char, initValue: bool) -> bool = false;
    addPriority(handle: RRHandle, eid: *const c_char, priority: *const c_char) -> bool = false;
    addPriorityNoRegen(handle: RRHandle, eid: *const c_char, priority: *const c_char) -> bool = false;
    addDelay(handle: RRHandle, eid: *const c_char, delay: *const c_char) -> bool = false;
    addDelayNoRegen(handle: RRHandle, eid: *const c_char, delay: *const c_char) -> bool = false;
    addEventAssignment(handle: RRHandle, eid: *const c_char, vid: *const c_char, formula: *const c_char) -> bool = false;
    addEventAssignmentNoRegen(handle: RRHandle, eid: *const c_char, vid: *const c_char, formula: *const c_char) -> bool = false;
    removeEventAssignments(handle: RRHandle, eid: *const c_char, vid: *const c_char) -> bool = false;
    removeEventAssignmentsNoRegen(handle: RRHandle, eid: *const c_char, vid: *const c_char) -> bool = false;
    removeEvent(handle: RRHandle, eid: *const c_char) -> bool = false;
    removeEventsNoRegen(handle: RRHandle, eid: *const c_char) -> bool = false;
    regenerate(handle: RRHandle) -> bool = false;
    validateCurrentSBML(handle: RRHandle) -> bool = false;
    // SBML utility methods
    getParamPromotedSBML(handle: RRHandle, sArg: *const c_char) -> *mut c_char = null_mut();
    setConfigurationXML(handle: RRHandle, caps: *const c_char) -> bool = false;
    getConfigurationXML(handle: RRHandle) -> *mut c_char = null_mut();
    getNumRegisteredIntegrators() -> c_int = -1;
    getRegisteredIntegratorName(n: c_int) -> *mut c_char = null_mut();
    getRegisteredIntegratorHint(n: c_int) -> *mut c_char = null_mut();
    getRegisteredIntegratorDescription(n: c_int) -> *mut c_char = null_mut();
    getNumInstantiatedIntegrators(handle: RRHandle) -> c_int = -1;
    setCurrentIntegrator(handle: RRHandle, nameOfIntegrator: *const c_char) -> c_int = -1;
    getCurrentIntegratorName(handle: RRHandle) -> *mut c_char = null_mut();
    getCurrentIntegratorDescription(handle: RRHandle) -> *mut c_char = null_mut();
    getCurrentIntegratorHint(handle: RRHandle) -> *mut c_char = null_mut();
    getNumberOfCurrentIntegratorParameters(handle: RRHandle) -> c_int = -1;
    getCurrentIntegratorNthParameterName(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentIntegratorNthParameterDescription(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentIntegratorNthParameterDisplayName(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentIntegratorNthParameterHint(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentIntegratorNthParameterType(handle: RRHandle, n: c_int) -> c_int = -1;
    resetCurrentIntegratorParameters(handle: RRHandle) -> c_int = -1;
    getListOfCurrentIntegratorParameterNames(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getCurrentIntegratorParameterDescription(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    getCurrentIntegratorParameterHint(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    getCurrentIntegratorParameterType(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    getCurrentIntegratorParameterInt(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    setCurrentIntegratorParameterInt(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int = -1;
    getCurrentIntegratorParameterUInt(handle: RRHandle, parameterName: *const c_char) -> c_uint = u32::MAX;
    setCurrentIntegratorParameterUInt(handle: RRHandle, parameterName: *const c_char, value: c_uint) -> c_int = -1;
    getCurrentIntegratorParameterDouble(handle: RRHandle, parameterName: *const c_char) -> f64 = -1.0;
    setCurrentIntegratorParameterDouble(handle: RRHandle, parameterName: *const c_char, value: f64) -> c_int = -1;
    getCurrentIntegratorParameterString(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    setCurrentIntegratorParameterString(handle: RRHandle, parameterName: *const c_char, value: *const c_char) -> c_int = -1;
    getCurrentIntegratorParameterBoolean(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    setCurrentIntegratorParameterBoolean(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int = -1;
    getCurrentIntegratorParameterDoubleArray(handle: RRHandle, parameterName: *const c_char, value: *mut *mut f64, len: *mut c_int) -> c_int = -1;
    setCurrentIntegratorParameterDoubleArray(handle: RRHandle, parameterName: *const c_char, value: *mut f64, len: c_int) -> c_int = -1;
    setCurrentIntegratorIndividualTolerance(handle: RRHandle, sid: *const c_char, value: f64) -> c_int = -1;
    getNumRegisteredSteadyStateSolvers() -> c_int = -1;
    getRegisteredSteadyStateSolverName(n: c_int) -> *mut c_char = null_mut();
    getRegisteredSteadyStateSolverHint(n: c_int) -> *mut c_char = null_mut();
    getRegisteredSteadyStateSolverDescription(n: c_int) -> *mut c_char = null_mut();
    setCurrentSteadyStateSolver(handle: RRHandle, nameOfSteadyStateSolver: *const c_char) -> c_int = -1;
    getCurrentSteadyStateSolverName(handle: RRHandle) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverDescription(handle: RRHandle) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverHint(handle: RRHandle) -> *mut c_char = null_mut();
    getNumberOfCurrentSteadyStateSolverParameters(handle: RRHandle) -> c_int = -1;
    getCurrentSteadyStateSolverNthParameterName(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverNthParameterDisplayName(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverNthParameterDescription(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverNthParameterHint(handle: RRHandle, n: c_int) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverNthParameterType(handle: RRHandle, n: c_int) -> c_int = -1;
    resetCurrentSteadyStateSolverParameters(handle: RRHandle) -> c_int = -1;
    getListOfCurrentSteadyStateSolverParameterNames(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getCurrentSteadyStateSolverParameterDescription(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverParameterHint(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    getCurrentSteadyStateSolverParameterType(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    getCurrentSteadyStateSolverParameterInt(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    setCurrentSteadyStateSolverParameterInt(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int = -1;
    getCurrentSteadyStateSolverParameterUInt(handle: RRHandle, parameterName: *const c_char) -> c_uint = u32::MAX;
    setCurrentSteadyStateSolverParameterUInt(handle: RRHandle, parameterName: *const c_char, value: c_uint) -> c_int = -1;
    getCurrentSteadyStateSolverParameterDouble(handle: RRHandle, parameterName: *const c_char) -> f64 = -1.0;
    setCurrentSteadyStateSolverParameterDouble(handle: RRHandle, parameterName: *const c_char, value: f64) -> c_int = -1;
    getCurrentSteadyStateSolverParameterString(handle: RRHandle, parameterName: *const c_char) -> *mut c_char = null_mut();
    setCurrentSteadyStateSolverParameterString(handle: RRHandle, parameterName: *const c_char, value: *const c_char) -> c_int = -1;
    getCurrentSteadyStateSolverParameterBoolean(handle: RRHandle, parameterName: *const c_char) -> c_int = -1;
    setCurrentSteadyStateSolverParameterBoolean(handle: RRHandle, parameterName: *const c_char, value: c_int) -> c_int = -1;
    setTimeStart(handle: RRHandle, timeStart: f64) -> bool = false;
    setTimeEnd(handle: RRHandle, timeEnd: f64) -> bool = false;
    setNumPoints(handle: RRHandle, numberOfPoints: c_int) -> bool = false;
    setTimes(handle: RRHandle, times: *const f64, size: c_int) -> bool = false;
    setTimeCourseSelectionList(handle: RRHandle, list: *const c_char) -> bool = false;
    getTimeCourseSelectionList(handle: RRHandle) -> *mut RRStringArray = null_mut();
    simulate(handle: RRHandle) -> *mut RRCData = null_mut();
    simulateNoReturn(handle: RRHandle) -> bool = false;
    getSimulationResult(handle: RRHandle) -> *mut RRCData = null_mut();
    getSimulationResultAsDoubleMatrix(handle: RRHandle) -> RRHandle = null_mut();
    simulateEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> *mut RRCData = null_mut();
    simulateExNoReturn(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> bool = false;
    simulateTimes(handle: RRHandle, times: *const f64, size: c_int) -> *mut RRCData = null_mut();
    oneStep(handle: RRHandle, currentTime: f64, stepSize: f64, value: *mut f64) -> bool = false;
    getTimeStart(handle: RRHandle, timeStart: *mut f64) -> bool = false;
    getTimeEnd(handle: RRHandle, timeEnd: *mut f64) -> bool = false;
    getNumPoints(handle: RRHandle, numPoints: *mut c_int) -> bool = false;
    steadyState(handle: RRHandle, value: *mut f64) -> bool = false;
    computeSteadyStateValues(handle: RRHandle) -> *mut RRVector = null_mut();
    setSteadyStateSelectionList(handle: RRHandle, list: *const c_char) -> bool = false;
    getSteadyStateSelectionList(handle: RRHandle) -> *mut RRStringArray = null_mut();
    // Get and Set Routines
    getValue(handle: RRHandle, symbolId: *const c_char, value: *mut f64) -> bool = false;
    setValue(handle: RRHandle, symbolId: *const c_char, value: f64) -> bool = false;
    getFloatingSpeciesConcentrations(handle: RRHandle) -> *mut RRVector = null_mut();
    getRateRuleValues(handle: RRHandle) -> *mut RRVector = null_mut();
    getFloatingSpeciesAmounts(handle: RRHandle) -> *mut RRVector = null_mut();
    getIndependentFloatingSpeciesConcentrations(handle: RRHandle) -> *mut RRVector = null_mut();
    getDependentFloatingSpeciesConcentrations(handle: RRHandle) -> *mut RRVector = null_mut();
    getIndependentFloatingSpeciesAmounts(handle: RRHandle) -> *mut RRVector = null_mut();
    getDependentFloatingSpeciesAmounts(handle: RRHandle) -> *mut RRVector = null_mut();
    getBoundarySpeciesConcentrations(handle: RRHandle) -> *mut RRVector = null_mut();
    getBoundarySpeciesAmounts(handle: RRHandle) -> *mut RRVector = null_mut();
    // Parameter Group
    getGlobalParameterValues(handle: RRHandle) -> *mut RRVector = null_mut();
    setBoundarySpeciesByIndex(handle: RRHandle, index: c_uint, value: f64) -> bool = false;
    setFloatingSpeciesInitialConcentrationByIndex(handle: RRHandle, index: c_int, value: f64) -> bool = false;
    getFloatingSpeciesInitialConcentrationByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    setFloatingSpeciesByIndex(handle: RRHandle, index: c_int, value: f64) -> bool = false;
    setGlobalParameterByIndex(handle: RRHandle, index: c_int, value: f64) -> bool = false;
    getBoundarySpeciesByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    getFloatingSpeciesByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    getGlobalParameterByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    getCompartmentByIndex(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    setCompartmentByIndex(handle: RRHandle, index: c_int, value: f64) -> bool = false;
    setFloatingSpeciesConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool = false;
    setBoundarySpeciesConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool = false;
    getFullJacobian(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getReducedJacobian(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getEigenvalues(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    // Stoichiometry methods
    getStoichiometryMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getLinkMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getNrMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getConservationMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    // Initial condition Methods
    setFloatingSpeciesInitialConcentrations(handle: RRHandle, vec: *mut RRVector) -> bool = false;
    getFloatingSpeciesInitialConcentrations(handle: RRHandle) -> *mut RRVector = null_mut();
    getFloatingSpeciesInitialConditionIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getFloatingSpeciesInitialConcentrationIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    // Reaction rates
    getNumberOfReactions(handle: RRHandle) -> c_int = -1;
    getReactionRate(handle: RRHandle, index: c_int, rate: *mut f64) -> bool = false;
    getReactionRates(handle: RRHandle) -> *mut RRVector = null_mut();
    getReactionRatesEx(handle: RRHandle, vec: *mut RRVector) -> *mut RRVector = null_mut();
    getRatesOfChange(handle: RRHandle) -> *mut RRVector = null_mut();
    getRatesOfChangeIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getRateOfChange(handle: RRHandle, index: c_int, value: *mut f64) -> bool = false;
    getRatesOfChangeEx(handle: RRHandle, vec: *mut RRVector) -> *mut RRVector = null_mut();
    evalModel(handle: RRHandle) -> bool = false;
    // Get number family
    getNumberOfCompartments(handle: RRHandle) -> c_int = -1;
    getNumberOfBoundarySpecies(handle: RRHandle) -> c_int = -1;
    getNumberOfFloatingSpecies(handle: RRHandle) -> c_int = -1;
    getNumberOfGlobalParameters(handle: RRHandle) -> c_int = -1;
    getNumberOfDependentSpecies(handle: RRHandle) -> c_int = -1;
    getNumberOfIndependentSpecies(handle: RRHandle) -> c_int = -1;
    // Get Ids family
    getReactionIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getBoundarySpeciesIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getBoundarySpeciesConcentrationIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getFloatingSpeciesIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getDependentFloatingSpeciesIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getIndependentFloatingSpeciesIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getFloatingSpeciesConcentrationIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getGlobalParameterIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getCompartmentIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getEigenvalueIds(handle: RRHandle) -> *mut RRStringArray = null_mut();
    getAvailableTimeCourseSymbols(handle: RRHandle) -> *mut RRList = null_mut();
    getAvailableSteadyStateSymbols(handle: RRHandle) -> *mut RRList = null_mut();
    // MCA methods
    getElasticityCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getUnscaledElasticityCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getUnscaledFluxControlCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getFluxControlCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getUnscaledConcentrationControlCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getConcentrationControlCoefficientIds(handle: RRHandle) -> *mut RRList = null_mut();
    getUnscaledElasticityMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getScaledElasticityMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getScaledFloatingSpeciesElasticity(handle: RRHandle, reactionId: *const c_char, speciesId: *const c_char, value: *mut f64) -> bool = false;
    getUnscaledConcentrationControlCoefficientMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getScaledConcentrationControlCoefficientMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getUnscaledFluxControlCoefficientMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getScaledFluxControlCoefficientMatrix(handle: RRHandle) -> *mut RRDoubleMatrix = null_mut();
    getuCC(handle: RRHandle, variable: *const c_char, parameter: *const c_char, value: *mut f64) -> bool = false;
    getCC(handle: RRHandle, variable: *const c_char, parameter: *const c_char, value: *mut f64) -> bool = false;
    getDiffStepSize(handle: RRHandle, value: *mut f64) -> bool = false;
    setDiffStepSize(handle: RRHandle, value: f64) -> bool = false;
    getEE(handle: RRHandle, name: *const c_char, species: *const c_char, value: *mut f64) -> bool = false;
    getuEE(handle: RRHandle, name: *const c_char, species: *const c_char, value: *mut f64) -> bool = false;
    // Stochastic methods
    getSeed(handle: RRHandle, seed: *mut i64) -> bool = false;
    setSeed(handle: RRHandle, seed: i64) -> bool = false;
    gillespie(handle: RRHandle) -> *mut RRCData = null_mut();
    gillespieEx(handle: RRHandle, timeStart: f64, timeEnd: f64) -> *mut RRCData = null_mut();
    gillespieOnGrid(handle: RRHandle) -> *mut RRCData = null_mut();
    gillespieOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int) -> *mut RRCData = null_mut();
    gillespieMeanOnGrid(handle: RRHandle, numberOfSimulations: c_int) -> *mut RRCData = null_mut();
    gillespieMeanOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfPoints: c_int, numberOfSimulations: c_int) -> *mut RRCData = null_mut();
    gillespieMeanSDOnGrid(handle: RRHandle, numberOfSimulations: c_int) -> *mut RRCData = null_mut();
    gillespieMeanSDOnGridEx(handle: RRHandle, timeStart: f64, timeEnd: f64, numberOfSteps: c_int, numberOfSimulations: c_int) -> *mut RRCData = null_mut();
    // Reset methods
    reset(handle: RRHandle) -> bool = false;
    resetAll(handle: RRHandle) -> bool = false;
    resetToOrigin(handle: RRHandle) -> bool = false;
    resetParameter(handle: RRHandle) -> bool = false;
    // RRPlugins
    _getNumIndFloatingSpecies(handle: RRHandle) -> c_int = -1;
    _getNumRateRules(handle: RRHandle) -> c_int = -1;
    _getTime(handle: RRHandle) -> f64 = -1.0;
    _getStateVector(handle: RRHandle) -> c_int = -1;
    _getStateVectorRate(handle: RRHandle, time: f64, value: *mut f64) -> () = ();
    setTimeCourseSelectionListEx(handle: RRHandle, length: c_int, list: *mut *const c_char) -> bool = false;
    setSteadyStateSelectionListEx(handle: RRHandle, length: c_int, list: *mut *const c_char) -> bool = false;
}
