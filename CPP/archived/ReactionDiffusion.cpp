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

	std::cout << "Matrix A: " << std::endl;
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

	// Allocate memory for f1 and f2
	f1 = new double[Nx*Ny];
	f2 = new double[Nx*Ny];

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
	std::cout << "\nTime to set initial conditions: " << InitialConditions.format() << " seconds" << std::endl;

	// // Print the initial conditions
	// std::cout << "U matrix (Initial):\n" << std::endl;
	// printFullMatrix(U, Nx, Ny);

	// std::cout << "V matrix (Intial):\n" << std::endl;
	// printFullMatrix(V, Nx, Ny);

}


// Reaction terms for the PDEs
// void ReactionDiffusion::solve_f1(){
// 	// return eps * u * (1.0 - u) * (u - (v + b)/a);
// 	for (int row = 0; row < Nx; ++row){
// 		for (int col = 0; col < Ny; ++col){
// 			f1[row*Ny + col] = eps * U[row*Ny + col] * (1.0 - U[row*Ny + col]) * (U[row*Ny + col] - (V[row*Ny + col] + b)/a);
// 		}
// 	}
// }

double ReactionDiffusion::solve_f1(double& u, double& v){
	// return eps * u * (1.0 - u) * (u - (v + b)/a);
	return eps * u * (1.0 - u) * (u - (v + b)/a);
}


// void ReactionDiffusion::solve_f2(){
// 	// return u * u * u - v;
// 	for (int row = 0; row < Nx; ++row){
// 		for (int col = 0; col < Ny; ++col){
// 			f2[row*Ny + col] = U[row*Ny + col] * U[row*Ny + col] * U[row*Ny + col] - V[row*Ny + col];
// 		}
// 	}
// }

double ReactionDiffusion::solve_f2(double& u, double& v){
	// return u * u * u - v;
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

	temp1 = new double[2*2];
	temp2 = new double[2*2];
	temp3 = new double[2*2];
	temp4 = new double[2*2];

	temp12 = new double[2*2];
	temp34 = new double[2*2];

	boost::timer::cpu_timer Solver;

	for (double t = dt; t <= T; t += dt){	
		
		// // Recalculate f1 and f2
		// solve_f1();
		// solve_f2();

		// Reset the dU and dV arrays
		for (int i = 0; i < Nx*Ny; i++){
			dU[i] = 0.0;
			dV[i] = 0.0;
		}

		// Each thread will have a private temp1 temp2 temp3 temp4 arrays of size 2x2 - these will be used to store the results of the twobytwo matrix multiplication
		// The result will be added to the dU and dV arrays after all the cells have been solved and the threads are rejoined

		// Iterate through the cells - note will need a loop over the rows and columns of the grid to Nx-1 and row*Ny-1
		// #pragma omp parallel for schedule(dynamic)
		for (int row = 0; row < Ny-1; ++row){
			// #pragma omp parallel for
			for (int cell = row*Nx; cell < (row+1)*Nx-1; ++cell){

				// For U
				TwoByTwo(A, 2, &U[cell], Nx, temp1);
				TwoByTwo(&U[cell], Nx, A, 2, temp2);

				// temp12 = temp1 + temp2
				for (int i = 0; i < 2*2; i++){
					temp12[i] = temp1[i] + temp2[i];
				}
				// F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, A, 2, &U[cell], Nx, 1.0, &dU[cell], Nx);
				// F77NAME(dgemm)("N", "N", 2, 2, 2, mu1, &U[cell], Nx, A, 2, 1.0,	&dU[cell], Nx);
				
				// For V
				TwoByTwo(A, 2, &V[cell], Nx, temp3);
				TwoByTwo(&V[cell], Nx, A, 2, temp4);

				// temp34 = temp3 + temp4
				for (int i = 0; i < 2*2; i++){
					temp34[i] = temp3[i] + temp4[i];
				}
				// F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, A, 2, &V[cell], Nx, 1.0, &dV[cell], Nx);
				// F77NAME(dgemm)("N", "N", 2, 2, 2, mu2, &V[cell], Nx, A, 2, 1.0, &dV[cell], Nx);

				// Add the results to the dU and dV arrays
				// #pragma omp parallel for

				# pragma omp critical
				for (int y_shift = 0; y_shift < 2; ++y_shift){
					for (int x_shift = 0; x_shift < 2; ++x_shift){
						dU[cell + x_shift + y_shift*Nx] += mu1 * temp12[x_shift + y_shift*2] + solve_f1(U[cell + x_shift + y_shift*Nx], V[cell + x_shift + y_shift*Nx]);
						dV[cell + x_shift + y_shift*Nx] += mu2 * temp34[x_shift + y_shift*2] + solve_f2(U[cell + x_shift + y_shift*Nx], V[cell + x_shift + y_shift*Nx]);
					}
				}

				// for (int i = 0; i < 2; ++i){
				// 	for (int j = 0; j < 2; ++j){
				// 		std::cout << "dU[" << cell + i*Ny + j << "] = " << dU[cell + i*Ny + j] << std::endl;
				// 	}
				// }
			
			// std::cout << "Row: " << row << std::endl;

			}

		}

		// Update the U and V arrays
		for (int node = 0; node < Nx*Ny; ++node){
			U[node] += dt * dU[node];
			V[node] += dt * dV[node];
		}

		std::cout << "Solver Time: " << t << std::endl;
		
	}
	
	printFullMatrix(U, Nx, Ny);
	std::cout << "Time to solve: " << Solver.format() << " seconds" << std::endl;
}






void ReactionDiffusion::writeToFile(){
	std::ofstream outfile;
	outfile.open("output.txt");

	// Print the solution to the file - format: x y u v
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			outfile << row*dx << " " << col*dy << " " << U[row*Ny + col] << " " << V[row*Ny + col] << std::endl;
		}
	}

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


// TwoByTwoMatrixMultiplication - multiplies two 2x2 matrices and returns the result in a 2x2 matrix - all matrices are stored in row major order
void TwoByTwo(double* A, int lda, double* B, int ldb, double* C){

	// Multiply the two matrices
	for (int col = 0; col < 2; ++col){
		for (int row = 0; row < 2; ++row){
			C[0] = A[0] * B[0] + A[1] * B[ldb];
			C[1] = A[0] * B[1] + A[1] * B[1 + ldb];
			C[2] = A[lda] * B[0] + A[1 + lda] * B[ldb];
			C[3] = A[lda] * B[1] + A[1 + lda] * B[1 + ldb];
		}
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
