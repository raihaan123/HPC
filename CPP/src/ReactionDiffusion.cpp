#include "ReactionDiffusion.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>

#include <omp.h>
// #include <cblas.h>

#define F77NAME(x) x##_

// Matrix matrix multiplication with BLAS
extern "C" {
	void F77NAME(dgemm)(char *transa, char *transb,
						const int& m, const int& n, const int& k,
						const double& alpha, double *a,
						const int& lda, double *b, const int& ldb,
						const double& beta, double *c, const int& ldc);
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

	// Allocate memory for A (2x2 general matrix)
	A = new double[2*2];

	// Building the fundamental shift matrix A
	A[0] = -1.0;
	A[1] = 1.0;
	A[2] = 1.0;
	A[3] = -1.0;

	printFullMatrix(A, 2, 2);

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

	// // Print the initial conditions
	// std::cout << "U matrix (Initial):\n" << std::endl;
	// printFullMatrix(U, Nx, Ny);

	// std::cout << "V matrix (Intial):\n" << std::endl;
	// printFullMatrix(V, Nx, Ny);

}


// Reaction terms for the PDEs
double ReactionDiffusion::f1(double& u, double& v){
	return eps * u * (1.0 - u) * (u - (v + b)/a);
}


double ReactionDiffusion::f2(double& u, double& v){
	return u * u * u - v;
}


/* Solver

--> The solution domain will be grouped into 'cells' of size dx*dy=1*1 - so 4 nodes per cell - so in total (Nx-1)*(Ny-1) cells
--> Each cell will be solved independently - these will be OpenMP tasks allocated to the team of threads
	--> f1 and f2 will be solved for each node in the cell - added to the dU and dV arrays
	--> A simplified Laplacian function will be applied to each node in the cell (AU + UA and AV + VA) - also added to the dU and dV arrays (reduced)
	-->	Once all the cells have been solved, the dU and dV arrays will update the U and V arrays respectively
--> The threads will be decoupled from any specific cell - they will be assigned to cells in a round-robin fashion

*/




void ReactionDiffusion::solve()
{

	boost::timer::cpu_timer Solver;

	for (int node = 0; node < Nx*Ny; ++node){
		// Calculate the reaction terms
		dU[node] = f1(U[node], V[node]);
		dV[node] = f2(U[node], V[node]);

		// Calculate the Laplacian terms - multiplication of the fundamental shift matrix with the U and V sub-matrices - allocated to several threads
	}

	// #pragma omp parallel
	// {
	// 	#pragma omp single
	// 	{
	// 		#pragma omp task
	// 		for (int cell = 0; cell < (Nx-1)*(Ny-1); ++cell){

	// 				// For U
	// 				F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, A, 2, &U[cell], Ny, 1.0, &dU[cell], Ny);
	// 			}

	// 		#pragma omp task
	// 		for (int cell = 0; cell < (Nx-1)*(Ny-1); ++cell){
				
	// 			// Reverse the order of matrix multiplication
	// 			F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, &U[cell], Ny, A, 2, 1.0, &dU[cell], Ny);
	// 		}

	// 		#pragma omp task
	// 		for (int cell = 0; cell < (Nx-1)*(Ny-1); ++cell){
	// 			// For V
	// 			F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, A, 2, &V[cell], Ny, 1.0, &dV[cell], Ny);
	// 		}

	// 		#pragma omp task
	// 		for (int cell = 0; cell < (Nx-1)*(Ny-1); ++cell){
	// 			// Reverse the order of matrix multiplication
	// 			F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, &V[cell], Ny, A, 2, 1.0, &dV[cell], Ny);
	// 		}
	// 	}
	// }

	// Non-parallel version
	#pragma omp parallel for schedule(static)
	for (int cell = 0; cell < (Nx-1)*(Ny-1); ++cell){

		// For U
		F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, A, 2, &U[cell], Ny, 1.0, &dU[cell], Ny);
		// Reverse the order of matrix multiplication
		F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, &U[cell], Ny, A, 2, 1.0, &dU[cell], Ny);

		// For V
		F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, A, 2, &V[cell], Ny, 1.0, &dV[cell], Ny);
		// Reverse the order of matrix multiplication
		F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, &V[cell], Ny, A, 2, 1.0, &dV[cell], Ny);
	}


	printFullMatrix(dV, Nx, Ny);
	std::cout << "Time to calculate reaction terms: " << Solver.format() << " seconds" << std::endl;





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



void GNUPlot(){
	/* Need to run the following commands in the terminal:
		set pm3d at st
		set view map
		set cbrange[0:1]
		splot 'output.txt' using 1:2:3 w l palette
	*/
	
	// Run the commands to the terminal 
	system("gnome-terminal -e 'bash -c \"set pm3d at st; set view map; set cbrange[0:1]; splot \'output.txt\' using 1:2:3 w l palette; set term x11; set output; pause -1;\"'");

}