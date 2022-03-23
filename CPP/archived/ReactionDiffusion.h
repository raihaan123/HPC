#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

class ReactionDiffusion{
	private:
		// Methods
		double solve_f1(double& u, double& v);
		double solve_f2(double& u, double& v);
		// void solve_f1();
		// void solve_f2();
		// double Laplacian(double* node);
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

		// f1 and f2
		double *f1, *f2;

		// Shift matrix A
		double *A;

		// Solver temp variables
		double* temp1;
		double* temp2;
		double* temp3;
		double* temp4;

		double* temp12;
		double* temp34;
		
	public:
		// Constructor
		ReactionDiffusion();

		// Destructor
		~ReactionDiffusion();

		void setParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps);		
		void setInitialConditions();
		void solve();
		void writeToFile();
		// void openmp_fun();
};


// Helper functions
void TwoByTwo(double* A, int lda, double* B, int ldb, double* C);
void printFullMatrix(double *A, int Nx, int Ny);
void printBandedSymmetricMatrix(double* A, int N);
void GNUPlot();

#endif