#include "ReactionDiffusion.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>


#include <omp.h>


// blas definitions
#define F77NAME(x) x##_

// Symmetric banded matrix A multiplication with a full matrix U - use dsbmv
extern "C" {

	void F77NAME(dsbmv)(const char *uplo, const int *n, const int *k,
						const double *alpha, const double *A, const int *lda,
						const double *x, const int *incx, const double *beta,
						double *y, const int *incy);
}

// Constructor with initializer list
ReactionDiffusion::ReactionDiffusion(){
	std::cout << "Constructing ReactionDiffusion object..." << std::endl;

}

// Destructor
ReactionDiffusion::~ReactionDiffusion(){
	std::cout << "Destructing ReactionDiffusion object..." << std::endl;

	// Free memory
	delete[] U;
	delete[] V;
	delete[] A;
}


void ReactionDiffusion::setParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps){
	// Set parameters to class attributes - 'this' is a pointer to the initialised object!
	this->dt = dt;
	this->T = T;
	this->Nx = Nx;
	this->Ny = Ny;
	this->a = a;
	this->b = b;
	this->mu1 = mu1;
	this->mu2 = mu2;
	this->eps = eps;

	// Allocate memory for U and V
	U = new double[Nx*Ny];
	V = new double[Nx*Ny];

	// Allocate memory for dU and dV
	dU = new double[Nx*Ny];
	dV = new double[Nx*Ny];

	// A is a square symmetric banded matrix of size Ny x Ny - store in banded form - has -2 on the diagonal and -1 on the off-diagonal
	double *A = new double[3*Ny];

	// Boost timer start
	boost::timer::cpu_timer FillingMatrix;

	// Populate the matrix A
    for (int i = 0; i < Ny-1; i++){
        A[i*3] = -2.0;
		A[i*3 + 1] = -1.0;
    }

	// Neumann boundary conditions
	A[0] = -1;
	A[(Ny-1)*3] = -1.0;

	// End timer
	std::cout << "Time to set A: " << FillingMatrix.format() << " seconds" << std::endl;

    // Print the A matrix
    for (int i = 0; i < Ny; i++){
        for (int j = 0; j < 3; j++){
            std::cout << A[i*3 + j] << " ";
        }
        std::cout << std::endl;
    }

	// Set initial conditions
	setInitialConditions();
}


// Set the initial conditions
void ReactionDiffusion::setInitialConditions(){
	// Set the initial conditions to U and V
	// For y>ly/2, U = 1, 0 everywhere else
	// For x<lx/2, V = a/2, 0 everywhere else
	// note x and y are the 'coordinates' of the grid points - not the actual x and y values - of U and V elements
	for (int i = 0; i < Nx*Ny; i++){
		int x = i % Nx;
		int y = i / Nx;

		if (y > Ny/2){
			U[i] = 1.0;
			V[i] = 0.0;
		}
		else{
			U[i] = 0.0;
			V[i] = 0.0;
		}

		if (x < Nx/2){
			U[i] = 0.0;
			V[i] = a/2.0;
		}
	}

	// Print the initial conditions
	std::cout << "Initial conditions:" << std::endl;
	// first u 
	for (int i = 0; i < Nx; i++){
		for (int j = 0; j < Ny; j++){
			std::cout << U[i*Ny + j] << " ";
		}
		std::cout << std::endl;
	}

}



// Solver
void ReactionDiffusion::solve(){
	// Print hi
	std::cout << "Hi!" << std::endl;
}


void ReactionDiffusion::writeToFile(){
	
	std::ofstream outfile;
	outfile.open("output.txt");

	// Print the solution to the file
	for (int i = 0; i < Nx*Ny; i++){
		outfile << U[i] << " ";
	}
	outfile << std::endl;

	// Close the file
	outfile.close();
}