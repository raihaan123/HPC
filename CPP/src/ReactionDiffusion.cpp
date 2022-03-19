#include "ReactionDiffusion.h"
#include <omp.h>
#include <iostream>

using namespace std;

// Constructor with initializer list
ReactionDiffusion::ReactionDiffusion(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps){
	double *U = new double[Nx*Ny];
	double *V = new double[Nx*Ny];

	// A is a square symmetric banded matrix of size Ny x Ny - store in banded form - has -2 on the diagonal and -1 on the off-diagonal
	double *A = new double[3*Ny];

    for (int i = 0; i < Ny; i++){
        A[i*3] = -2.0;
    }

    // Set second column = -1
    for (int i = 0; i < Ny; i++){
        A[i*3 + 1] = -1.0;
    }

    // Print the A matrix
    for (int i = 0; i < Ny; i++){
        for (int j = 0; j < 3; j++){
            std::cout << A[i*3 + j] << " ";
        }
        std::cout << std::endl;
    }


}

// Solver
void ReactionDiffusion::solve(){
}