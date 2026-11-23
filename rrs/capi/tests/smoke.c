/* Smoke test for the roadrunner_nitro C API (libroadrunner rrc_api.h compatible).
 * Build & run from the repo root, after `cargo build --release -p roadrunner_nitro_c_api` in rrs/:
 *   cc -Wall -I rrs/capi/include rrs/capi/tests/smoke.c -L rrs/target/release -lroadrunner_nitro_c_api \
 *      -Wl,-rpath,rrs/target/release -o /tmp/rrc_smoke && /tmp/rrc_smoke
 * Leak check (macOS): leaks --atExit -- /tmp/rrc_smoke   (expects "0 leaks")
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
    EQS(getlibSBMLVersion(NULL), "none");
    EQS(getAPIVersion(), "1.0");
    assert(strcmp(solverTypeToString(8), "double") == 0);
    char* loc = getRRCAPILocation();
    assert(loc && strlen(loc) > 0);
    freeText(loc);

    /* NULL handles fail cleanly */
    assert(!loadSBML(NULL, "<sbml/>"));
    assert(strstr(getLastError(), "Failed to cast to a valid RoadRunner handle"));
    assert(simulate(NULL) == NULL && getNumberOfReactions(NULL) == -1);

    /* a model: S1 -> S2 (k1), S2 -> S3 (k2), $X0 -> S1 (v0) */
    static const char* sbml =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\"><model id=\"m\">"
        "<listOfCompartments><compartment id=\"c\" size=\"1\" constant=\"true\"/></listOfCompartments>"
        "<listOfSpecies>"
        "<species id=\"X0\" compartment=\"c\" initialConcentration=\"1\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"true\" constant=\"false\"/>"
        "<species id=\"S1\" compartment=\"c\" initialConcentration=\"0\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\" constant=\"false\"/>"
        "<species id=\"S2\" compartment=\"c\" initialConcentration=\"0\" hasOnlySubstanceUnits=\"false\" boundaryCondition=\"false\" constant=\"false\"/>"
        "</listOfSpecies>"
        "<listOfParameters><parameter id=\"v0\" value=\"2\" constant=\"true\"/><parameter id=\"k1\" value=\"0.5\" constant=\"true\"/>"
        "<parameter id=\"k2\" value=\"0.25\" constant=\"true\"/></listOfParameters>"
        "<listOfReactions>"
        "<reaction id=\"J0\" reversible=\"false\"><listOfReactants><speciesReference species=\"X0\" stoichiometry=\"1\" constant=\"true\"/></listOfReactants>"
        "<listOfProducts><speciesReference species=\"S1\" stoichiometry=\"1\" constant=\"true\"/></listOfProducts>"
        "<kineticLaw><math xmlns=\"http://www.w3.org/1998/Math/MathML\"><apply><times/><ci>v0</ci><ci>X0</ci></apply></math></kineticLaw></reaction>"
        "<reaction id=\"J1\" reversible=\"false\"><listOfReactants><speciesReference species=\"S1\" stoichiometry=\"1\" constant=\"true\"/></listOfReactants>"
        "<listOfProducts><speciesReference species=\"S2\" stoichiometry=\"1\" constant=\"true\"/></listOfProducts>"
        "<kineticLaw><math xmlns=\"http://www.w3.org/1998/Math/MathML\"><apply><times/><ci>k1</ci><ci>S1</ci></apply></math></kineticLaw></reaction>"
        "<reaction id=\"J2\" reversible=\"false\"><listOfReactants><speciesReference species=\"S2\" stoichiometry=\"1\" constant=\"true\"/></listOfReactants>"
        "<kineticLaw><math xmlns=\"http://www.w3.org/1998/Math/MathML\"><apply><times/><ci>k2</ci><ci>S2</ci></apply></math></kineticLaw></reaction>"
        "</listOfReactions></model></sbml>";
    RRHandle rr = createRRInstance();
    assert(rr);
    assert(setTempFolder(rr, "/tmp/rrc"));
    EQS(getTempFolder(rr), "/tmp/rrc");
    assert(!isModelLoaded(rr) && getNumberOfReactions(rr) == -1 && hasError());
    assert(loadSBML(rr, sbml) && isModelLoaded(rr));
    assert(getNumberOfReactions(rr) == 3 && getNumberOfFloatingSpecies(rr) == 2);
    assert(setTimeCourseSelectionList(rr, "time, [S1], S2"));
    RRCDataPtr res = simulateEx(rr, 0, 10, 11);
    assert(res && res->RSize == 11 && res->CSize == 3);
    assert(strcmp(res->ColumnHeaders[1], "[S1]") == 0);
    assert(getRRCDataElement(res, 10, 0, &d) && d == 10);
    /* S1(t) = 4 (1 - exp(-t/2)) */
    assert(getRRCDataElement(res, 10, 1, &d) && fabs(d - 4 * (1 - exp(-5.0))) < 1e-5);
    freeRRCData(res);
    assert(steadyState(rr, &d) && d < 1e-6);
    assert(getValue(rr, "[S1]", &d) && fabs(d - 4) < 1e-6);   /* v0 / k1 */
    assert(getValue(rr, "S2", &d) && fabs(d - 8) < 1e-6);     /* v0 / k2 */
    /* d S1 / d v0 = S1 / v0: scaled 1 */
    assert(getCC(rr, "S1", "v0", &d) && fabs(d - 1) < 1e-6);
    assert(getCC(rr, "J2", "k1", &d) && fabs(d) < 1e-6);
    assert(!getCC(rr, "S1", "nope", &d) && hasError());
    RRDoubleMatrixPtr jac = getFullJacobian(rr);
    assert(jac && jac->RSize == 2 && jac->Data[0] == -0.5);
    freeMatrix(jac);
    RRVectorPtr fs = getFloatingSpeciesConcentrations(rr);
    assert(fs && fs->Count == 2);
    freeVector(fs);
    RRStringArrayPtr ids = getFloatingSpeciesIds(rr);
    EQS(stringArrayToString(ids), "S1 S2");
    freeStringArray(ids);
    RRListPtr syms = getAvailableTimeCourseSymbols(rr);
    assert(getListLength(syms) == 12);
    freeRRList(syms);
    assert(freeRRInstance(rr));
    assert(!freeRRInstance(NULL));

    printf("smoke OK\n");
    return 0;
}
