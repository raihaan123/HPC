#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

class ReactionDiffusion{
	private:
		// Methods
		void TimeIntegrate();

		// Parameters
		double dt;
		double T;
		int Nx;
		int Ny;
		double a;
		double b;
		double mu1;
		double mu2;
		double eps;

		// dx and dy are const
		const double dx = 1.0;
		const double dy = 1.0;


		// U and V are double matrices of size Nx x Ny - represent the solution of the PDEs
		double *U, *V;

		// Derivative of U and V
		double *dU, *dV;

		// Shift matrix A
		double *A;

		
	public:
		// Constructor
		ReactionDiffusion();

		// Destructor
		~ReactionDiffusion();

		void setParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps);		
		void setInitialConditions();
		void solve();
		void writeToFile();
};

#endif