#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

class ReactionDiffusion{
	private:
		// U and V are double matrices of size Nx x Ny - represent the solution of the PDEs
		double *U, *V;

		// Derivative of U and V
		double *dU, *dV;

		// Shift matrix A
		double *A;

		
	public:
		// Constructor
		ReactionDiffusion(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps);

		void TimeIntegrate();
		void SetParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps);
		void SetInitialConditions(double *U, double *V);

		void solve();
};

#endif