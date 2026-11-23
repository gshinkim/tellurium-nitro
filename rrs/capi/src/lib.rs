//! libroadrunner C API (wrappers/C: rrc_api.h, rrc_types.h, rrc_utilities.h) backed by roadrunner_nitro.
#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals, clippy::missing_safety_doc)]

pub mod api;
pub mod types;
pub mod utilities;

use std::cell::RefCell;
use std::ffi::CString;

thread_local! {
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
}

/// Record the last error for this thread (read back by `hasError` / `getLastError`).
pub fn set_error(msg: impl Into<String>) {
    let msg = CString::new(msg.into().replace('\0', " ")).unwrap();
    LAST_ERROR.with(|e| *e.borrow_mut() = Some(msg));
}

pub(crate) fn has_error() -> bool {
    LAST_ERROR.with(|e| e.borrow().is_some())
}

/// Pointer into the thread-local store; valid until the next `set_error` on this thread.
pub(crate) fn last_error_ptr() -> *mut libc::c_char {
    LAST_ERROR.with(|e| e.borrow().as_ref().map_or(c"No Error".as_ptr(), |s| s.as_ptr()) as *mut _)
}
