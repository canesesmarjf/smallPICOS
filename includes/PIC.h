#ifndef H_PIC
#define H_PIC

// System libraries:
#include <iostream>
#include <cmath>
#include <vector>

/*
#ifdef __linux__
#include <chrono>//C++11 standard
#include <random>//C++11 standard

#elif __APPLE__
//define something else for apple compiler
#endif
*/

// User-defined libraries:
#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"

// Parallelization libraries:
// #include <omp.h>
// #include "mpi_main.h"

using namespace std;
using namespace arma;

class PIC_TYP
{

protected:

	// MPI methods:
	/*
	void MPI_AllreduceVec(const params_TYP * params, arma::vec * v);

	void MPI_SendVec(const params_TYP * params, arma::vec * v);

	void MPI_ReduceVec(const params_TYP * params, arma::vec * v);

	void MPI_Allgathervec(const params_TYP * params, arma::vec * field);

	void MPI_Recvvec(const params_TYP * params, arma::vec * field);

	void MPI_Recv_AllFields(const params_TYP * params, fields_TYP * fields);
	*/

	// Ghost contributions:
	void fillGhosts(arma::vec * C);

	void fill4Ghosts(arma::vec * v);

	void fillGhost_AllFields(const params_TYP * params, fields_TYP * fields);

	// Smoothing:
	void smooth(arma::vec * v, double as);

	// PIC related:
	void interpolateFields(const params_TYP * params, ions_TYP * IONS, const fields_TYP * fields);

	void interpolateElectrons(const params_TYP * params, ions_TYP * IONS, const electrons_TYP * electrons);

	void interpolateScalarField(const params_TYP * params, ions_TYP * IONS, const arma::vec * F_m, arma::vec * F_p);

	/*
	void interpEM(const params_TYP * params, const mesh_TYP * mesh, const ions_TYP * IONS, const fields_TYP * fields, arma::rowvec * ZN, arma::rowvec * EM);

	void calculateF(const params_TYP * params, const ions_TYP * IONS, arma::rowvec * ZN, arma::rowvec * EM, arma::rowvec * F);
	*/

	void eim(const params_TYP * params, fields_TYP * fields, ions_TYP * IONS);

	void calculateIonMoments(const params_TYP * params, fields_TYP * fields, ions_TYP * IONS);

	void calculateDerivedIonMoments(const params_TYP * params, ions_TYP * IONS);

  public:

	PIC_TYP(const params_TYP * params, const mesh_TYP * mesh,fields_TYP * fields, vector<ions_TYP> * IONS, electrons_TYP * electrons);

	void assignCell(const params_TYP * params, const mesh_TYP * mesh, ions_TYP * IONS);

	/*
	void advanceParticles(const params_TYP * params, const mesh_TYP * mesh, fields_TYP * fields, vector<ions_TYP> * IONS);
	*/
	
	void assignCell_AllSpecies(const params_TYP * params, const mesh_TYP * mesh, vector<ions_TYP> * IONS);

	void interpolateFields_AllSpecies(const params_TYP * params, vector<ions_TYP> * IONS, const fields_TYP * fields);

	void interpolateElectrons_AllSpecies(const params_TYP * params, vector<ions_TYP> * IONS, const electrons_TYP * electrons);

	void extrapolateMoments_AllSpecies(const params_TYP * params, fields_TYP * fields, vector<ions_TYP> * IONS);

};

#endif
