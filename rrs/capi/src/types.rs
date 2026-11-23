//! `#[repr(C)]` mirrors of rrc_types.h plus allocation helpers.
//!
//! Everything handed to C is allocated with libc `malloc` and released with libc `free`
//! (by the `free*` functions in utilities.rs), so C callers may also `free()` strings directly.
use libc::{c_char, c_int, c_void};
use std::ffi::CStr;
use std::mem::size_of;
use std::ptr::null_mut;

pub type RRHandle = *mut c_void;

#[repr(C)]
pub struct RRVector {
    pub Count: c_int,
    pub Data: *mut f64,
}

#[repr(C)]
pub struct RRStringArray {
    pub Count: c_int,
    pub String: *mut *mut c_char,
}

#[repr(C)]
pub struct RRDoubleMatrix {
    pub RSize: c_int,
    pub CSize: c_int,
    pub Data: *mut f64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RRComplex {
    pub re: f64,
    pub imag: f64,
}

#[repr(C)]
pub struct RRComplexVector {
    pub Count: c_int,
    pub Data: *mut RRComplex,
}

#[repr(C)]
pub struct RRComplexMatrix {
    pub RSize: c_int,
    pub CSize: c_int,
    pub Data: *mut RRComplex,
}

#[repr(C)]
pub struct RRCData {
    pub RSize: c_int,
    pub CSize: c_int,
    pub Data: *mut f64,
    pub Weights: *mut f64,
    pub ColumnHeaders: *mut *mut c_char,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ListItemType {
    litString,
    litInteger,
    litDouble,
    litList,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RRParameterType {
    ptString,
    ptBool,
    ptInteger,
    ptDouble,
    ptVector,
    ptMatrix,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub union RRListItemData {
    pub iValue: c_int,
    pub dValue: f64,
    pub sValue: *mut c_char,
    pub lValue: *mut RRList,
}

#[repr(C)]
pub struct RRListItem {
    pub ItemType: ListItemType,
    pub data: RRListItemData,
}

#[repr(C)]
pub struct RRList {
    pub Count: c_int,
    pub Items: *mut *mut RRListItem,
}

#[cfg(target_pointer_width = "64")]
const _: () = assert!(
    size_of::<RRVector>() == 16 && size_of::<RRStringArray>() == 16 && size_of::<RRDoubleMatrix>() == 16
        && size_of::<RRComplex>() == 16 && size_of::<RRComplexVector>() == 16 && size_of::<RRComplexMatrix>() == 16
        && size_of::<RRCData>() == 32 && size_of::<RRListItem>() == 16 && size_of::<RRList>() == 16
        && size_of::<ListItemType>() == 4
);

// ---------------------------------------------------------------- allocation helpers

/// malloc one `T` and move `v` into it.
pub fn alloc<T>(v: T) -> *mut T {
    let p = unsafe { libc::malloc(size_of::<T>()) } as *mut T;
    assert!(!p.is_null(), "malloc failed");
    unsafe { p.write(v) };
    p
}

/// malloc a copy of `s` (NULL when empty).
pub fn alloc_slice<T: Copy>(s: &[T]) -> *mut T {
    if s.is_empty() {
        return null_mut();
    }
    let p = unsafe { libc::malloc(std::mem::size_of_val(s)) } as *mut T;
    assert!(!p.is_null(), "malloc failed");
    unsafe { p.copy_from_nonoverlapping(s.as_ptr(), s.len()) };
    p
}

/// malloc'd NUL-terminated copy of `s` (free with `freeText` / `free`).
pub fn text(s: &str) -> *mut c_char {
    let p = unsafe { libc::malloc(s.len() + 1) } as *mut c_char;
    assert!(!p.is_null(), "malloc failed");
    unsafe {
        p.copy_from_nonoverlapping(s.as_ptr() as *const c_char, s.len());
        *p.add(s.len()) = 0;
    }
    p
}

/// Borrow a C string argument (None for NULL; invalid UTF-8 is replaced).
pub unsafe fn cstr(p: *const c_char) -> Option<String> {
    (!p.is_null()).then(|| CStr::from_ptr(p).to_string_lossy().into_owned())
}

fn strings<S: AsRef<str>>(items: &[S]) -> *mut *mut c_char {
    alloc_slice(&items.iter().map(|s| text(s.as_ref())).collect::<Vec<_>>())
}

pub fn vector(data: &[f64]) -> *mut RRVector {
    alloc(RRVector { Count: data.len() as c_int, Data: alloc_slice(data) })
}

/// Row-major `rows x cols` matrix.
pub fn matrix(rows: usize, cols: usize, data: &[f64]) -> *mut RRDoubleMatrix {
    assert_eq!(rows * cols, data.len());
    alloc(RRDoubleMatrix { RSize: rows as c_int, CSize: cols as c_int, Data: alloc_slice(data) })
}

pub fn string_array<S: AsRef<str>>(items: &[S]) -> *mut RRStringArray {
    alloc(RRStringArray { Count: items.len() as c_int, String: strings(items) })
}

/// Row-major simulation result; `headers.len()` must equal the column count.
pub fn cdata<S: AsRef<str>>(rows: usize, data: &[f64], weights: Option<&[f64]>, headers: &[S]) -> *mut RRCData {
    assert_eq!(rows * headers.len(), data.len());
    alloc(RRCData {
        RSize: rows as c_int,
        CSize: headers.len() as c_int,
        Data: alloc_slice(data),
        Weights: weights.map_or(null_mut(), alloc_slice),
        ColumnHeaders: strings(headers),
    })
}
