//! rrc_utilities.h: text, list, vector/matrix helpers, free routines and *ToString (upstream formatting).
use crate::set_error;
use crate::types::*;
use libc::{c_char, c_int, free};
use std::ffi::CStr;
use std::fmt::Write;
use std::ptr::null_mut;

/// C++ `ostream << double` default formatting (`%g`, precision 6).
pub(crate) fn g(x: f64) -> String {
    let mut buf = [0 as c_char; 64];
    unsafe {
        libc::snprintf(buf.as_mut_ptr(), buf.len(), c"%g".as_ptr(), x);
        CStr::from_ptr(buf.as_ptr()).to_string_lossy().into_owned()
    }
}

unsafe fn s(p: *const c_char) -> String {
    cstr(p).unwrap_or_default()
}

unsafe fn slice<'a, T>(p: *const T, n: c_int) -> &'a [T] {
    if p.is_null() || n <= 0 { &[] } else { std::slice::from_raw_parts(p, n as usize) }
}

/// Bounds-checked row-major index into a `rows x cols` block; sets the upstream error text on failure.
fn index(what: &str, rows: c_int, cols: c_int, r: c_int, c: c_int) -> Option<usize> {
    if r < 0 || c < 0 || r >= rows || c >= cols {
        set_error(format!("Index out range in {what}: {r}, {c}"));
        return None;
    }
    Some((r * cols + c) as usize)
}

// ---------------------------------------------------------------- text

#[no_mangle]
pub unsafe extern "C" fn createText(t: *const c_char) -> *mut c_char {
    if t.is_null() { null_mut() } else { libc::strdup(t) }
}

#[no_mangle]
pub extern "C" fn createTextMemory(count: c_int) -> *mut c_char {
    if count < 0 { null_mut() } else { unsafe { libc::calloc(count as usize + 1, 1) as *mut c_char } }
}

#[no_mangle]
pub unsafe extern "C" fn getFileContent(fName: *const c_char) -> *mut c_char {
    let name = s(fName);
    match std::fs::read_to_string(&name) {
        Ok(c) => text(&c),
        Err(e) => {
            set_error(format!("RoadRunner exception: could not read file '{name}': {e}"));
            null_mut()
        }
    }
}

// ---------------------------------------------------------------- lists

#[no_mangle]
pub extern "C" fn createRRList() -> *mut RRList {
    alloc(RRList { Count: 0, Items: null_mut() })
}

#[no_mangle]
pub unsafe extern "C" fn freeRRList(list: *mut RRList) {
    let Some(l) = list.as_ref() else { return };
    for &item in slice(l.Items, l.Count) {
        if let Some(it) = item.as_ref() {
            match it.ItemType {
                ListItemType::litString => free(it.data.sValue as _),
                ListItemType::litList => freeRRList(it.data.lValue),
                _ => {}
            }
        }
        free(item as _);
    }
    free(l.Items as _);
    free(list as _);
}

#[no_mangle]
pub unsafe extern "C" fn getListLength(myList: *mut RRList) -> c_int {
    myList.as_ref().map_or_else(|| { set_error("List argument is null in getListLength"); -1 }, |l| l.Count)
}

#[no_mangle]
pub extern "C" fn createIntegerItem(value: c_int) -> *mut RRListItem {
    alloc(RRListItem { ItemType: ListItemType::litInteger, data: RRListItemData { iValue: value } })
}

#[no_mangle]
pub extern "C" fn createDoubleItem(value: f64) -> *mut RRListItem {
    alloc(RRListItem { ItemType: ListItemType::litDouble, data: RRListItemData { dValue: value } })
}

#[no_mangle]
pub unsafe extern "C" fn createStringItem(value: *const c_char) -> *mut RRListItem {
    alloc(RRListItem { ItemType: ListItemType::litString, data: RRListItemData { sValue: createText(value) } })
}

#[no_mangle]
pub extern "C" fn createListItem(value: *mut RRList) -> *mut RRListItem {
    alloc(RRListItem { ItemType: ListItemType::litList, data: RRListItemData { lValue: value } })
}

/// Appends `*item` to `list`; the list takes ownership. Returns its index, or -1.
#[no_mangle]
pub unsafe extern "C" fn addItem(list: *mut RRList, item: *mut *mut RRListItem) -> c_int {
    let (Some(l), false) = (list.as_mut(), item.is_null()) else {
        set_error("Null argument in addItem()");
        return -1;
    };
    let n = l.Count.max(0) as usize;
    let items = libc::realloc(l.Items as _, (n + 1) * size_of::<*mut RRListItem>()) as *mut *mut RRListItem;
    if items.is_null() {
        set_error("Failed allocating memory in addItem()");
        return -1;
    }
    *items.add(n) = *item;
    l.Items = items;
    l.Count = n as c_int + 1;
    n as c_int
}

#[no_mangle]
pub unsafe extern "C" fn getListItem(list: *mut RRList, index: c_int) -> *mut RRListItem {
    list.as_ref().and_then(|l| slice(l.Items, l.Count).get(usize::try_from(index).ok()?).copied()).unwrap_or(null_mut())
}

unsafe fn is(item: *mut RRListItem, t: ListItemType) -> bool {
    item.as_ref().is_some_and(|i| i.ItemType == t)
}

#[no_mangle]
pub unsafe extern "C" fn isListItemInteger(item: *mut RRListItem) -> c_int {
    is(item, ListItemType::litInteger) as c_int
}

#[no_mangle]
pub unsafe extern "C" fn isListItemDouble(item: *mut RRListItem) -> c_int {
    is(item, ListItemType::litDouble) as c_int
}

#[no_mangle]
pub unsafe extern "C" fn isListItemString(item: *mut RRListItem) -> c_int {
    is(item, ListItemType::litString) as c_int
}

#[no_mangle]
pub unsafe extern "C" fn isListItemList(item: *mut RRListItem) -> c_int {
    is(item, ListItemType::litList) as c_int
}

#[no_mangle]
pub unsafe extern "C" fn isListItem(item: *mut RRListItem, itemType: ListItemType) -> c_int {
    is(item, itemType) as c_int
}

#[no_mangle]
pub unsafe extern "C" fn getIntegerListItem(item: *mut RRListItem, value: *mut c_int) -> c_int {
    if !is(item, ListItemType::litInteger) || value.is_null() {
        return 0;
    }
    *value = (*item).data.iValue;
    1
}

#[no_mangle]
pub unsafe extern "C" fn getDoubleListItem(item: *mut RRListItem, value: *mut f64) -> c_int {
    if !is(item, ListItemType::litDouble) || value.is_null() {
        return 0;
    }
    *value = (*item).data.dValue;
    1
}

#[no_mangle]
pub unsafe extern "C" fn getStringListItem(item: *mut RRListItem) -> *mut c_char {
    if is(item, ListItemType::litString) { (*item).data.sValue } else { null_mut() }
}

#[no_mangle]
pub unsafe extern "C" fn getList(item: *mut RRListItem) -> *mut RRList {
    if is(item, ListItemType::litList) { (*item).data.lValue } else { null_mut() }
}

// ---------------------------------------------------------------- free routines

#[no_mangle]
pub unsafe extern "C" fn freeRRCData(handle: *mut RRCData) -> c_int {
    if let Some(d) = handle.as_ref() {
        free(d.Data as _);
        free(d.Weights as _);
        if !d.ColumnHeaders.is_null() {
            for &h in slice(d.ColumnHeaders, d.CSize) {
                free(h as _);
            }
        }
        free(d.ColumnHeaders as _);
        free(handle as _);
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn freeText(text: *mut c_char) -> c_int {
    free(text as _);
    1
}

#[no_mangle]
pub unsafe extern "C" fn freeStringArray(sl: *mut RRStringArray) -> c_int {
    if let Some(a) = sl.as_ref() {
        for &p in slice(a.String, a.Count) {
            free(p as _);
        }
        free(a.String as _);
        free(sl as _);
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn freeVector(vector: *mut RRVector) -> c_int {
    if let Some(v) = vector.as_ref() {
        free(v.Data as _);
        free(vector as _);
    }
    1
}

#[no_mangle]
pub unsafe extern "C" fn freeMatrix(matrix: *mut RRDoubleMatrix) -> c_int {
    if let Some(m) = matrix.as_ref() {
        free(m.Data as _);
        free(matrix as _);
    }
    1
}

// ---------------------------------------------------------------- vectors & matrices

#[no_mangle]
pub unsafe extern "C" fn getVectorLength(vector: *mut RRVector) -> c_int {
    vector.as_ref().map_or_else(|| { set_error("Vector argument is null in getVectorLength"); -1 }, |v| v.Count)
}

/// Zero-initialised; Data is NULL for size 0.
#[no_mangle]
pub extern "C" fn createVector(size: c_int) -> *mut RRVector {
    if size < 0 {
        set_error(format!("Negative size in createVector: {size}"));
        return null_mut();
    }
    vector(&vec![0.0; size as usize])
}

#[no_mangle]
pub unsafe extern "C" fn getVectorElement(vector: *mut RRVector, index: c_int, value: *mut f64) -> c_int {
    let Some(v) = vector.as_ref() else {
        set_error("Vector argument is null in getVectorElement");
        return 0;
    };
    if index < 0 || index >= v.Count || value.is_null() {
        set_error(format!("Index out range in getVectorElement: {index}"));
        return 0;
    }
    *value = *v.Data.add(index as usize);
    1
}

#[no_mangle]
pub unsafe extern "C" fn setVectorElement(vector: *mut RRVector, index: c_int, value: f64) -> c_int {
    let Some(v) = vector.as_ref() else {
        set_error("Vector argument is null in setVectorElement");
        return 0;
    };
    if index < 0 || index >= v.Count {
        set_error(format!("Index out range in setVectorElement: {index}"));
        return 0;
    }
    *v.Data.add(index as usize) = value;
    1
}

/// Zero-initialised r x c matrix; NULL (with error) when r*c == 0.
#[no_mangle]
pub extern "C" fn createRRMatrix(r: c_int, c: c_int) -> *mut RRDoubleMatrix {
    if r <= 0 || c <= 0 {
        set_error("Dimensions for new RRDoubleMatrix in createRRMatrix are zero");
        return null_mut();
    }
    matrix(r as usize, c as usize, &vec![0.0; r as usize * c as usize])
}

#[no_mangle]
pub unsafe extern "C" fn getMatrixNumRows(m: *mut RRDoubleMatrix) -> c_int {
    m.as_ref().map_or_else(|| { set_error("Matrix argument is null in getMatrixNumRows"); -1 }, |m| m.RSize)
}

#[no_mangle]
pub unsafe extern "C" fn getMatrixNumCols(m: *mut RRDoubleMatrix) -> c_int {
    m.as_ref().map_or_else(|| { set_error("Matrix argument is null in getMatrixNumCols"); -1 }, |m| m.CSize)
}

#[no_mangle]
pub unsafe extern "C" fn getMatrixElement(m: *mut RRDoubleMatrix, r: c_int, c: c_int, value: *mut f64) -> c_int {
    let Some(m) = m.as_ref() else {
        set_error("Matrix argument is null in getMatrixElement");
        return 0;
    };
    let Some(i) = index("getMatrixElement", m.RSize, m.CSize, r, c) else { return 0 };
    if value.is_null() {
        return 0;
    }
    *value = *m.Data.add(i);
    1
}

#[no_mangle]
pub unsafe extern "C" fn setMatrixElement(m: *mut RRDoubleMatrix, r: c_int, c: c_int, value: f64) -> c_int {
    let Some(m) = m.as_ref() else {
        set_error("Matrix argument is null in setMatrixElement");
        return 0;
    };
    let Some(i) = index("setMatrixElement", m.RSize, m.CSize, r, c) else { return 0 };
    *m.Data.add(i) = value;
    1
}

#[no_mangle]
pub unsafe extern "C" fn getComplexMatrixElement(m: *mut RRComplexMatrix, r: c_int, c: c_int, value: *mut RRComplex) -> c_int {
    let Some(m) = m.as_ref() else {
        set_error("Matrix argument is null in getComplexMatrixElement");
        return 0;
    };
    let Some(i) = index("getComplexMatrixElement", m.RSize, m.CSize, r, c) else { return 0 };
    if value.is_null() {
        return 0;
    }
    *value = *m.Data.add(i);
    1
}

#[no_mangle]
pub unsafe extern "C" fn setComplexMatrixElement(m: *mut RRComplexMatrix, r: c_int, c: c_int, value: *mut RRComplex) -> c_int {
    let Some(m) = m.as_ref() else {
        set_error("Matrix argument is null in setComplexMatrixElement");
        return 0;
    };
    let Some(i) = index("setComplexMatrixElement", m.RSize, m.CSize, r, c) else { return 0 };
    if value.is_null() {
        return 0;
    }
    *m.Data.add(i) = *value;
    1
}

// ---------------------------------------------------------------- RRCData

#[no_mangle]
pub unsafe extern "C" fn getRRDataNumRows(rrData: *mut RRCData) -> c_int {
    rrData.as_ref().map_or_else(|| { set_error("result argument is null in getResultNumRows"); -1 }, |d| d.RSize)
}

#[no_mangle]
pub unsafe extern "C" fn getRRDataNumCols(rrData: *mut RRCData) -> c_int {
    rrData.as_ref().map_or_else(|| { set_error("result argument is null in getResultNumCols"); -1 }, |d| d.CSize)
}

#[no_mangle]
pub unsafe extern "C" fn getRRCDataElement(rrData: *mut RRCData, r: c_int, c: c_int, value: *mut f64) -> c_int {
    let Some(d) = rrData.as_ref() else {
        set_error("result argument is null in getResultElement");
        return 0;
    };
    let Some(i) = index("getResultElement", d.RSize, d.CSize, r, c) else { return 0 };
    if value.is_null() {
        return 0;
    }
    *value = *d.Data.add(i);
    1
}

/// Returns the internal label pointer (do not free).
#[no_mangle]
pub unsafe extern "C" fn getRRDataColumnLabel(rrData: *mut RRCData, column: c_int) -> *mut c_char {
    let Some(d) = rrData.as_ref() else {
        set_error("result argument is null in getResultColumnLabel");
        return null_mut();
    };
    if column < 0 || column >= d.CSize || d.ColumnHeaders.is_null() {
        set_error(format!("Index out range in getResultColumnLabel: {column}"));
        return null_mut();
    }
    *d.ColumnHeaders.add(column as usize)
}

// ---------------------------------------------------------------- to-string (byte-compatible with upstream)

#[no_mangle]
pub unsafe extern "C" fn vectorToString(vec: *mut RRVector) -> *mut c_char {
    let Some(v) = vec.as_ref() else {
        set_error("Null vector in vectorToString");
        return null_mut();
    };
    let mut out: String = slice(v.Data, v.Count).iter().map(|&x| g(x) + "\t").collect();
    out.push('\n');
    text(&out)
}

#[no_mangle]
pub unsafe extern "C" fn complexVectorToString(vec: *mut RRComplexVector) -> *mut c_char {
    let Some(v) = vec.as_ref() else {
        set_error("Null vector in vectorToString");
        return null_mut();
    };
    let mut out: String = slice(v.Data, v.Count).iter().map(|z| format!("({},{})\t", g(z.re), g(z.imag))).collect();
    out.push('\n');
    text(&out)
}

#[no_mangle]
pub unsafe extern "C" fn rrCDataToString(rrData: *mut RRCData) -> *mut c_char {
    let Some(d) = rrData.as_ref() else { return null_mut() };
    let mut out = (0..d.CSize as usize)
        .map(|i| if d.ColumnHeaders.is_null() { String::new() } else { s(*d.ColumnHeaders.add(i)) })
        .collect::<Vec<_>>()
        .join("\t");
    out.push('\n');
    let data = slice(d.Data, d.RSize * d.CSize);
    for row in data.chunks(d.CSize.max(1) as usize) {
        out += &row.iter().map(|&x| g(x)).collect::<Vec<_>>().join("\t");
        out.push('\n');
    }
    text(&out)
}

#[no_mangle]
pub unsafe extern "C" fn matrixToString(mat: *mut RRDoubleMatrix) -> *mut c_char {
    let Some(m) = mat.as_ref() else { return null_mut() };
    let mut out = String::from("\n");
    for row in slice(m.Data, m.RSize * m.CSize).chunks(m.CSize.max(1) as usize) {
        row.iter().for_each(|&x| { let _ = write!(out, "{}\t", g(x)); });
        out.push('\n');
    }
    text(&out)
}

#[no_mangle]
pub unsafe extern "C" fn complexMatrixToString(mat: *mut RRComplexMatrix) -> *mut c_char {
    let Some(m) = mat.as_ref() else { return null_mut() };
    let mut out = String::from("\n");
    for row in slice(m.Data, m.RSize * m.CSize).chunks(m.CSize.max(1) as usize) {
        row.iter().for_each(|z| { let _ = write!(out, "({},{})\t", g(z.re), g(z.imag)); });
        out.push('\n');
    }
    text(&out)
}

#[no_mangle]
pub unsafe extern "C" fn getNumberOfStringElements(list: *mut RRStringArray) -> c_int {
    list.as_ref().map_or(-1, |l| l.Count)
}

/// Returns a malloc'd copy (free with freeText).
#[no_mangle]
pub unsafe extern "C" fn getStringElement(list: *mut RRStringArray, index: c_int) -> *mut c_char {
    let Some(l) = list.as_ref() else { return null_mut() };
    if index < 0 || index >= l.Count {
        set_error("Index out of range");
        return null_mut();
    }
    createText(*l.String.add(index as usize))
}

#[no_mangle]
pub unsafe extern "C" fn stringArrayToString(list: *mut RRStringArray) -> *mut c_char {
    let Some(l) = list.as_ref() else { return null_mut() };
    text(&slice(l.String, l.Count).iter().map(|&p| s(p)).collect::<Vec<_>>().join(" "))
}

unsafe fn list_string(l: &RRList) -> String {
    let parts: Vec<String> = slice(l.Items, l.Count)
        .iter()
        .map(|&item| match item.as_ref() {
            None => String::new(),
            Some(it) => match it.ItemType {
                ListItemType::litString => format!("\"{}\"", s(it.data.sValue)),
                ListItemType::litInteger => it.data.iValue.to_string(),
                ListItemType::litDouble => g(it.data.dValue),
                ListItemType::litList => it.data.lValue.as_ref().map_or("{}".into(), |l| list_string(l)),
            },
        })
        .collect();
    format!("{{{}}}", parts.join(","))
}

#[no_mangle]
pub unsafe extern "C" fn listToString(list: *mut RRList) -> *mut c_char {
    list.as_ref().map_or(null_mut(), |l| text(&list_string(l)))
}

// ---------------------------------------------------------------- handle-based utilities

#[no_mangle]
pub unsafe extern "C" fn writeRRData(handle: RRHandle, fileNameAndPath: *const c_char) -> c_int {
    crate::api::write_rr_data(handle, fileNameAndPath)
}

/// Upstream is a no-op that returns true.
#[no_mangle]
pub extern "C" fn compileSource(_handle: RRHandle, _sourceFileNameAndPath: *const c_char) -> c_int {
    1
}
