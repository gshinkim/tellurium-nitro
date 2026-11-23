/* Smoke test for the roadrunner_nitro C API (libroadrunner rrc_api.h compatible).
 * Build & run from the repo root, after `cargo build --release -p roadrunner_nitro_c_api` in rrs/:
 *   cc -Wall -I rrs/capi/include rrs/capi/tests/smoke.c -L rrs/target/release -lroadrunner_nitro_c_api \
 *      -Wl,-rpath,rrs/target/release -o /tmp/rrc_smoke && /tmp/rrc_smoke
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rrc_api.h"

#define EQS(a, b) do { char* _s = (a); if (strcmp(_s, (b)) != 0) { printf("got [%s] want [%s]\n", _s, (b)); assert(0); } freeText(_s); } while (0)

int main(void)
{
    double d;
    int i;

    /* version + errors */
    char* v = getVersionStr();
    assert(v && strlen(v) > 0);
    printf("version %s (%d)\n", v, getVersion());
    freeText(v);
    assert(strcmp(getLastError(), "No Error") == 0 && !hasError());

    /* vectors: Rust allocates, C reads/writes the struct directly */
    RRVectorPtr vec = createVector(3);
    assert(getVectorLength(vec) == 3 && vec->Count == 3);
    vec->Data[0] = 1.5;
    assert(setVectorElement(vec, 1, 2.0) && setVectorElement(vec, 2, 1e-7));
    assert(getVectorElement(vec, 0, &d) && d == 1.5);
    assert(!getVectorElement(vec, 3, &d) && hasError());
    assert(strcmp(getLastError(), "Index out range in getVectorElement: 3") == 0);
    EQS(vectorToString(vec), "1.5\t2\t1e-07\t\n");
    freeVector(vec);
    assert(getVectorLength(NULL) == -1);

    /* matrices */
    RRDoubleMatrixPtr m = createRRMatrix(2, 3);
    assert(getMatrixNumRows(m) == 2 && getMatrixNumCols(m) == 3);
    m->Data[1 * m->CSize + 2] = 3.25;
    assert(getMatrixElement(m, 1, 2, &d) && d == 3.25);
    assert(setMatrixElement(m, 0, 0, 1.0/3.0));
    assert(!setMatrixElement(m, 2, 0, 1.0));
    EQS(matrixToString(m), "\n0.333333\t0\t0\t\n0\t0\t3.25\t\n");
    freeMatrix(m);
    assert(createRRMatrix(0, 3) == NULL);

    /* complex matrix (caller-owned) */
    struct RRComplex cells[2] = {{1, 2}, {0, 0}}, z = {5, -1}, out;
    struct RRComplexMatrix cm = {1, 2, cells};
    assert(setComplexMatrixElement(&cm, 0, 1, &z) && getComplexMatrixElement(&cm, 0, 1, &out));
    assert(out.re == 5 && out.imag == -1);
    EQS(complexMatrixToString(&cm), "\n(1,2)\t(5,-1)\t\n");

    /* lists: {123,{3.14159},"hi"} */
    RRListPtr l = createRRList();
    RRListItemPtr it = createIntegerItem(123);
    assert(addItem(l, &it) == 0);
    RRListPtr inner = createRRList();
    it = createDoubleItem(3.1415926);
    addItem(inner, &it);
    it = createListItem(inner);
    assert(addItem(l, &it) == 1);
    it = createStringItem("hi");
    assert(addItem(l, &it) == 2);
    assert(getListLength(l) == 3);
    assert(getListItem(l, 0)->data.iValue == 123 && isListItemInteger(getListItem(l, 0)));
    assert(getIntegerListItem(getListItem(l, 0), &i) && i == 123);
    assert(!getDoubleListItem(getListItem(l, 0), &d));
    assert(getDoubleListItem(getListItem(getList(getListItem(l, 1)), 0), &d) && d == 3.1415926);
    assert(strcmp(getStringListItem(getListItem(l, 2)), "hi") == 0 && isListItem(getListItem(l, 2), litString));
    assert(getListItem(l, 3) == NULL);
    EQS(listToString(l), "{123,{3.14159},\"hi\"}");
    freeRRList(l);

    /* string arrays (via the config-key store) */
    assert(setConfigBool("B_KEY", 1) && setConfigInt("A_KEY", 42) && setConfigDouble("C_KEY", 0.5));
    assert(getConfigBool("B_KEY") == 1 && getConfigInt("A_KEY") == 42 && getConfigDouble("C_KEY") == 0.5);
    assert(getConfigInt("NOPE") == -1);
    RRStringArrayPtr keys = getListOfConfigKeys();
    assert(getNumberOfStringElements(keys) == 3);
    EQS(getStringElement(keys, 1), "B_KEY");
    assert(getStringElement(keys, 3) == NULL);
    EQS(stringArrayToString(keys), "A_KEY B_KEY C_KEY");
    freeStringArray(keys);

    /* RRCData allocated by C with malloc, freed by the library (same allocator) */
    RRCDataPtr r = malloc(sizeof *r);
    r->RSize = 2; r->CSize = 2;
    r->Data = malloc(4 * sizeof(double));
    r->Weights = NULL;
    r->ColumnHeaders = malloc(2 * sizeof(char*));
    r->ColumnHeaders[0] = createText("time");
    r->ColumnHeaders[1] = createText("[S1]");
    double vals[4] = {0, 10, 0.5, 9.75};
    memcpy(r->Data, vals, sizeof vals);
    assert(getRRDataNumRows(r) == 2 && getRRDataNumCols(r) == 2);
    assert(getRRCDataElement(r, 1, 1, &d) && d == 9.75);
    assert(strcmp(getRRDataColumnLabel(r, 1), "[S1]") == 0);
    EQS(rrCDataToString(r), "time\t[S1]\n0\t10\n0.5\t9.75\n");
    freeRRCData(r);

    /* text + folders */
    char* t = createTextMemory(4);
    assert(t && t[0] == 0);
    freeText(t);
    assert(setTempFolder(NULL, "/tmp/rrc"));
    EQS(getTempFolder(NULL), "/tmp/rrc");
    EQS(getlibSBMLVersion(NULL), "none");
    EQS(getAPIVersion(), "1.0");
    assert(strcmp(solverTypeToString(8), "double") == 0);
    char* loc = getRRCAPILocation();
    assert(loc && strlen(loc) > 0);
    freeText(loc);

    /* stubs fail cleanly */
    assert(createRRInstance() == NULL);
    assert(!loadSBML(NULL, "<sbml/>"));
    assert(strstr(getLastError(), "not implemented yet"));
    assert(simulate(NULL) == NULL && getNumberOfReactions(NULL) == -1);

    printf("smoke OK\n");
    return 0;
}
