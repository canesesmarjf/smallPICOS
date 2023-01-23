#include<iostream>
#include<vector>

#define ARMA_ALLOW_FAKE_GCC

#include<armadillo>
#include<cmath>
#include <ctime>
#include <utility>

// User defined header files:
#include "types.h"
#include "initialize.h"

// Namespaces used in code:
using namespace std;
using namespace arma;

// Constants:
# define PI arma::datum::pi

int main(int argc, char* argv[])
{
    // Input parameters for simulation:
    params_TYP params;

    // Ion species vector:
    vector<ions_TYP> IONS;

    // Electron species object:
    electrons_TYP electrons;

    // Electromagnetic fields:
    fields_TYP fields;

    // Mesh object:
    mesh_TYP mesh;

    // Object to hold IC condition profiles:
    IC_TYP IC;

    // Initialization object:
    init_TYP init(&params, argc, argv);

    // =========================================================================
    // - Read "input_file.input" into "params"
    init.read_inputFile(&params);

    // - Read "ions_properties.ion" into "params":
    init.read_ionsPropertiesFile(&params);

    // - Create mesh using input parameters:
    init.create_mesh(&params, &mesh);

    // - Read IC profiles from external files into "IC":
    init.read_IC_profiles(&params, &mesh, &IC);

    // - Interpolate IC profiles to mesh grid:
    init.interpolate_IC_profiles(&params, &mesh, &IC);

    // - Initialize fields using IC field profiles:
    init.initialize_fields(&params,&IC,&fields);

    // - Initialize electrons using IC profiles:
    init.initialize_electrons(&params,&IC,&electrons);

    // Calculate particle weight initial condition profile:
    // PUT before initialize fields and electrons
    init.calculate_IC_particleWeight(&params,&IC,&mesh,&fields,&IONS);

    // fields.Bx_m.print("Bx_m");

    // - Initialize ions using IC profiles:
    //init.initialize_ions(&params,&IC,&mesh,&IONS);

    // Calculate global quantities: N_R, N_SP:
    //init.calculate_globalQuantities(&params,&mesh,&fields,&IONS);
}
