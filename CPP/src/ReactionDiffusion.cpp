#include "ReactionDiffusion.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>

#include <omp.h>
#include <cblas.h>


// blas definitions
#define F77NAME(x) x##_

// Symmetric banded matrix A multiplication with a full matrix U - use dsbmv!
extern "C" {
	void F77NAME(dsbmv)(const char *uplo, const int *n, const int *k,
						const double *alpha, const double *A, const int *lda,
						const double *x, const int *incx, const double *beta,
						double *y, const int *incy);

	void F77NAME(dgemm)(const char *transa, const char *transb, const int *m, const int *n, const int *k,
						const double *alpha, const double *A, const int *lda,
						const double *B, const int *ldb, const double *beta,
						double *C, const int *ldc);
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
	double *A = new double[2*Ny];

	// Full size A matrix
	double *A_full = new double[Ny*Ny];

	// Boost timer start
	boost::timer::cpu_timer FillingMatrix;

	/*
	Populate the matrix A!
	----------------------

	First row is all 1s from 2 to Ny
	Second row is -1 at start and finish, and -2 at the middle

	*/

	// Superdiagonal
	A[0] = 0.0;
	for (int col = 1; col < Ny; col++){
		A[col] = 1.0;
	}
	
	// Main diagonal
	for (int col = 0; col < Ny-1; col++){
		A[Ny + col] = -2.0;
	}

	// Neumann boundary conditions
	A[Ny] = -1.0;
	A[Ny*2 - 1] = -1.0;

	// End timer
	std::cout << "Time to set A: " << FillingMatrix.format() << " seconds" << std::endl;

    // Print the A matrix (nicely formatted!)
	// std::cout << "A matrix:" << std::endl;
    // printBandedSymmetricMatrix(A, Ny);

	// A_full is the full matrix A - main diagonal is all -2s and 1st super diagonal is all 1s - all else 0 - symmetric
	boost::timer::cpu_timer FillingMatrix2;

	for (int row = 0; row < Ny; row++){
		for (int col = 0; col < Ny; col++){
			if (row == col){
				A_full[row*Ny + col] = -2.0;
			}
			else if (row == col - 1 || row == col + 1){
				A_full[row*Ny + col] = 1.0;
			}
			else{
				A_full[row*Ny + col] = 0.0;
			}
		}
	}

	// Neumann boundary conditions
	A_full[0] = -1.0;
	A_full[Ny*Ny-1] = -1.0;

	// End timer
	std::cout << "Time to set A_full: " << FillingMatrix2.format() << " seconds" << std::endl;

	// Print the A_full matrix
	// std::cout << "A_full matrix:" << std::endl;
	// printFullMatrix(A_full, Ny, Ny);

}


// Set the initial conditions
void ReactionDiffusion::setInitialConditions(){
	/* 
	Set the initial conditions to U and V
	-------------------------------------

	--> For y>ly/2, U = 1, 0 everywhere else
	--> For x<lx/2, V = a/2, 0 everywhere else

	*/

	// Lx and Ly are always integers - Nx and Ny are odd integers
	const int Ly = (Ny-1)*dx;
	const int Lx = (Nx-1)*dx;

	const int row_bound = Ly/2;		// The first row of U equal to 0 (from the top)
	const int col_bound = Lx/2;		// The first column of V equal to 0 (from the left)

	// Allocate memory for U and V - all zero
	U = new double[Nx*Ny];
	V = new double[Nx*Ny];

	// All zero elements
	for (int i = 0; i < Nx*Ny; i++){
		U[i] = 0.0;
		V[i] = 0.0;
	}
	
	// Set the initial conditions
	boost::timer::cpu_timer InitialConditions;
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			if (row < row_bound){
				U[row*Ny + col] = 1.0;
			}

			if (col < col_bound){
				V[row*Ny + col] = a/2.0;
			}
		}
	}
	std::cout << "Time to set initial conditions: " << InitialConditions.format() << " seconds" << std::endl;

	// Print the initial conditions
	std::cout << "Initial conditions:" << std::endl;
	//print u
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			std::cout << U[row*Ny + col] << " ";
		}
		std::cout << std::endl;
	}
	// print v
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			std::cout << V[row*Ny + col] << " ";
		}
		std::cout << std::endl;
	}

}


// Reaction terms for the PDEs
void ReactionDiffusion::f1(){
	// epsilon * u * (1-u) * (u - (v+b)/a)
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			dU[row*Ny + col] += eps * U[row*Ny + col] * (1.0 - U[row*Ny + col]) * (U[row*Ny + col] - (V[row*Ny + col] + b)/a);
			std::cout <<"dU[" << row*Ny + col << "] = " << dU[row*Ny + col] << std::endl;
		}
	}
}

void ReactionDiffusion::f2(){
	// u**3 - v
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			dV[row*Ny + col] += U[row*Ny + col]*U[row*Ny + col]*U[row*Ny + col] - V[row*Ny + col];
			std::cout <<"dV[" << row*Ny + col << "] = " << dV[row*Ny + col] << std::endl;
		}
	}
}


// Solver
void ReactionDiffusion::solve(){
	std::cout << "Called solver!" << std::endl;

	// Allocate memory for dU and dV
	dU = new double[Nx*Ny];
	dV = new double[Nx*Ny];

	// for loop from 0 to T in steps of dt
	for (int t = 0; t <= T; t += dt){

		// Evaluate the reaction terms - methods with save results in dU and dV
		f1();
		f2();

		// mu1*(A*U + U*A) + f1
		// mu2*(A*V + V*A) + f2

		// k is the number of super-diagonals in the matrix A
		const int k = 1;

		// F77NAME(dsbmv)("U", &Ny, &k, &mu1, A, &Ny, U, &inc, &zero, dU, &inc);

		// Using cblas_dsbmv for one row of U at a time (not the whole matrix) - using a for loop
		for (int row = 0; row < Nx; ++row){

			// Note that beta = 1 because we are adding to dU (f1 is already in dU)
			std::cout << "Gulp!" << std::endl;
			cblas_dgemv(CblasRowMajor, CblasNoTrans, Ny, k, mu1, A_full, Ny, &U[row*Ny], 1, 1, &dU[row*Ny], 1);
			std::cout << "survived dsbmv" << std::endl;

			// same for v
			cblas_dgemv(CblasRowMajor, CblasNoTrans, Ny, k, mu2, A_full, Ny, &V[row*Ny], 1, 1, &dV[row*Ny], 1);
		
		}

		// Update U
		for (int row = 0; row < Nx; ++row){
			for (int col = 0; col < Ny; ++col){
				U[row*Ny + col] += dt * dU[row*Ny + col];
			}
		}
		std::cout << "survived U update" << std::endl;

		// Update V
		for (int row = 0; row < Nx; ++row){
			for (int col = 0; col < Ny; ++col){
				V[row*Ny + col] += dt * dV[row*Ny + col];
			}
		}
		std::cout << "survived V update" << std::endl;

	}
	

	//print the solution
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			std::cout << U[row*Ny + col] << " ";
		}
		std::cout << std::endl;
	}

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




/* Helper functions! */

void printFullMatrix(double* A, int Nx, int Ny){
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			std::cout << A[row*Ny + col] << " ";
		}
		std::cout << std::endl;
	}
}

void printBandedSymmetricMatrix(double* A, int N){
	for (int row = 0; row < 2; row++){
        for (int col = 0; col < N; col++){
            std::cout << A[row*N + col] << " ";
        }
        std::cout << std::endl;
    }

}