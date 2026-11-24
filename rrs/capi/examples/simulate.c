/* Simulate a model, find its steady state and a control coefficient through libroadrunner's C API,
 * on Tellurium Nitro's engine (libroadrunner_nitro_c_api). Build and run from the repo root:
 *
 *   cd rrs && cargo build --release -p roadrunner_nitro_c_api && cd ..
 *   cc -I rrs/capi/include rrs/capi/examples/simulate.c -L rrs/target/release -lroadrunner_nitro_c_api \
 *      -Wl,-rpath,$PWD/rrs/target/release -o simulate
 *   ./simulate                  # the built-in model
 *   ./simulate my_model.xml     # or any SBML file
 */
#include <stdio.h>
#include "rrc_api.h"

static const char *MODEL =
    "<?xml version='1.0' encoding='UTF-8'?>"
    "<sbml xmlns='http://www.sbml.org/sbml/level3/version2/core' level='3' version='2'><model id='m'>"
    "<listOfCompartments><compartment id='c' size='1' constant='true'/></listOfCompartments>"
    "<listOfSpecies>"
    "<species id='X0' compartment='c' initialConcentration='10' hasOnlySubstanceUnits='false' boundaryCondition='true' constant='false'/>"
    "<species id='S1' compartment='c' initialConcentration='0' hasOnlySubstanceUnits='false' boundaryCondition='false' constant='false'/>"
    "</listOfSpecies>"
    "<listOfParameters><parameter id='k1' value='0.5' constant='true'/><parameter id='k2' value='0.2' constant='true'/></listOfParameters>"
    "<listOfReactions>"
    "<reaction id='J1' reversible='false'><listOfReactants><speciesReference species='X0' stoichiometry='1' constant='true'/></listOfReactants>"
    "<listOfProducts><speciesReference species='S1' stoichiometry='1' constant='true'/></listOfProducts>"
    "<kineticLaw><math xmlns='http://www.w3.org/1998/Math/MathML'><apply><times/><ci>k1</ci><ci>X0</ci></apply></math></kineticLaw></reaction>"
    "<reaction id='J2' reversible='false'><listOfReactants><speciesReference species='S1' stoichiometry='1' constant='true'/></listOfReactants>"
    "<kineticLaw><math xmlns='http://www.w3.org/1998/Math/MathML'><apply><times/><ci>k2</ci><ci>S1</ci></apply></math></kineticLaw></reaction>"
    "</listOfReactions></model></sbml>";

int main(int argc, char **argv) {
    RRHandle rr = createRRInstance();
    if (!(argc > 1 ? loadSBMLFromFile(rr, argv[1]) : loadSBML(rr, MODEL))) {
        fprintf(stderr, "load failed: %s\n", getLastError());
        return 1;
    }

    /* time course: 0 to 10, 6 points */
    RRCDataPtr result = simulateEx(rr, 0, 10, 6);
    if (!result) { fprintf(stderr, "simulate failed: %s\n", getLastError()); return 1; }
    int rows = getRRDataNumRows(result), cols = getRRDataNumCols(result);
    for (int c = 0; c < cols; c++) printf("%12s", getRRDataColumnLabel(result, c));
    printf("\n");
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) { double v; getRRCDataElement(result, r, c, &v); printf("%12.6g", v); }
        printf("\n");
    }
    freeRRCData(result);

    /* steady state, then values and a control coefficient there (built-in model only) */
    double closeness;
    if (!steadyState(rr, &closeness)) { fprintf(stderr, "steady state failed: %s\n", getLastError()); return 1; }
    if (argc == 1) {
        double s1, cc;
        getValue(rr, "S1", &s1);
        getCC(rr, "J2", "k1", &cc);
        printf("steady state S1 = %g (expected 25), control coefficient of J2 on k1 = %g (expected 1)\n", s1, cc);
    }
    freeRRInstance(rr);
    return 0;
}
