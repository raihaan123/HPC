#include "ReactionDiffusion.h"
#include <iostream>
#include <fstream>
#include <boost/timer/timer.hpp>

#include <omp.h>


// Constructor
ReactionDiffusion::ReactionDiffusion(){
	std::cout << "\nConstructing ReactionDiffusion object...\n" << std::endl;
}


// Destructor
ReactionDiffusion::~ReactionDiffusion(){
	std::cout << "\nDestructing ReactionDiffusion object..." << std::endl;

	// Free memory
	delete[] U;
	delete[] V;
}


// Set parameters for the solver
void ReactionDiffusion::setParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps){
	// Set parameters to class attributes - 'this' is a pointer to the initialised object!
	this->dt = dt;
	this->T = T;
	this->Nx = Nx;
	this->Ny = Ny;
	this->a = 1/a;
	this->b = b;
	this->mu1 = mu1 * dt;       // mu1 is multiplied by dt to avoid repeated multiplication in the code
	this->mu2 = mu2 * dt;       // mu2 is multiplied by dt to avoid repeated multiplication in the code
	this->eps = eps;

	// Allocate memory for U and V
	U = new double[Nx*Ny];
	V = new double[Nx*Ny];

	// Allocate memory for dU and dV
	dU = new double[Nx*Ny];
	dV = new double[Nx*Ny];

    // Allocate memory for f1 and f2
    f1 = new double[Nx*Ny];
    f2 = new double[Nx*Ny];

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
	
	// Set the initial conditions
	boost::timer::cpu_timer InitialConditions;
	for (int row = 0; row < Nx; ++row){
		for (int col = 0; col < Ny; ++col){
			if (row < row_bound){
				U[row*Ny + col] = 1.0;
			}

			if (col < col_bound){
				V[row*Ny + col] = 1/(2*a);			// Note that a has been inverted previously
			}
		}
	}

	std::cout << "\nTime to set initial conditions: " << InitialConditions.format() << " seconds" << std::endl;

}


double ReactionDiffusion::solve_f1(double& u, double& v){
	// return eps * u * (1.0 - u) * (u - (v + b)/a);
	return dt * eps * u * (1.0 - u) * (u - (v + b)*a);
}


double ReactionDiffusion::solve_f2(double& u, double& v){
	// return u * u * u - v;
	return dt * (u * u * u - v);
}


void ReactionDiffusion::solve()
{
	boost::timer::cpu_timer Solver;

	for (double t = dt; t <= T; t += dt){

        // OpenMP parallelisation - each thread solves a different part of the grid - use tasks

        // Solve the diffusion equation - iterate over all interior nodes

        // Top left corner
        dU[0] = -2.0*U[0] + U[1] + U[Nx];
        dV[0] = -2.0*V[0] + V[1] + V[Nx];

        // Top edge
        // #pragma omp parallel for
        for (int col = 1; col < Nx-1; ++col){
            dU[col] = -3.0*U[col] + U[col-1] + U[col+1] + U[col+Nx];
            dV[col] = -3.0*V[col] + V[col-1] + V[col+1] + V[col+Nx];
        }

        // Top right corner
        dU[Nx-1] = -2.0*U[Nx-1] + U[Nx-2] + U[Nx*2-1];
        dV[Nx-1] = -2.0*V[Nx-1] + V[Nx-2] + V[Nx*2-1];


        #pragma omp parallel for
        for (int row = 1; row < Ny-1; ++row){
            #pragma omp parallel for
            for (int inner_node = row*Nx+1; inner_node < (row+1)*Nx-1; ++inner_node){

                // Apply the Laplacian stencil to the interior nodes
                dU[inner_node] = -4.0*U[inner_node] + U[inner_node-1] + U[inner_node+1] + U[inner_node-Nx] + U[inner_node+Nx];
                dV[inner_node] = -4.0*V[inner_node] + V[inner_node-1] + V[inner_node+1] + V[inner_node-Nx] + V[inner_node+Nx];
            }
        }

        // Left and right edges - memory is non-contiguous - U and V are row-major, so can solve both edges in parallel
        #pragma omp parallel for
        for (int row = 1; row < Ny-1; ++row){
            dU[row*Nx] = -3.0*U[row*Nx] + U[(row+1)*Nx] + U[(row-1)*Nx] + U[row*Nx+1];     // Left edge
            dV[row*Nx] = -3.0*V[row*Nx] + V[(row+1)*Nx] + V[(row-1)*Nx] + V[row*Nx+1];     // Left edge

            dU[(row+1)*Nx-1] = -3.0*U[(row+1)*Nx-1] + U[(row+1)*Nx-2] + U[(row)*Nx-1] + U[(row+2)*Nx-1];     // Right edge
            dV[(row+1)*Nx-1] = -3.0*V[(row+1)*Nx-1] + V[(row+1)*Nx-2] + V[(row)*Nx-1] + V[(row+2)*Nx-1];     // Right edge
        }


        // Bottom left corner
        dU[Nx*(Ny-1)] = -2.0*U[Nx*(Ny-1)] + U[Nx*(Ny-1)+1] + U[Nx*(Ny-2)];
        dV[Nx*(Ny-1)] = -2.0*V[Nx*(Ny-1)] + V[Nx*(Ny-1)+1] + V[Nx*(Ny-2)];

        // Bottom edge
        // #pragma omp parallel for
        for (int col = (Ny-1)*Nx+1; col < Ny*Nx-1; ++col){
            dU[col] = -3.0*U[col] + U[col-1] + U[col+1] + U[col-Ny];
            dV[col] = -3.0*V[col] + V[col-1] + V[col+1] + V[col-Ny];
        }

        // Bottom right corner
        dU[Nx*Ny-1] = -2.0*U[Nx*Ny-1] + U[Nx*Ny-2] + U[Nx*(Ny-1)-1];
        dV[Nx*Ny-1] = -2.0*V[Nx*Ny-1] + V[Nx*Ny-2] + V[Nx*(Ny-1)-1];


        // Update the U and V values
        #pragma omp parallel for 
        for (int node = 0; node < Nx*Ny; ++node){
            // U[node] += mu1*dU[node] + f1[node];
            // V[node] += mu2*dV[node] + f2[node];
			U[node] += mu1*dU[node] + solve_f1(U[node], V[node]);
			V[node] += mu2*dV[node] + solve_f2(U[node], V[node]);
        }

        // std::cout << "Simulation time: " << t << std::endl;
    }







	
	std::cout << "Time to solve: " << Solver.format() << " seconds" << std::endl;

}



void ReactionDiffusion::writeToFile(){
	std::ofstream outfile;
	outfile.open("output.txt");

	// Print the solution to the file - format: [x y u v]
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
