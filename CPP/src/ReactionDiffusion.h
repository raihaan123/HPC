#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#define dx 1.0
#define dy 1.0

/**
 * @class ReactionDiffusion
 * @brief Solve the reaction-diffusion equation using the finite difference method and in parallel
 */
class ReactionDiffusion{
	public:
		/// @brief Default constructor
		ReactionDiffusion();

		/// @brief Default destructor
		~ReactionDiffusion();
		
		/**
		 * @brief Set Parameters
		 * @details Set the parameters for the solver
		 * @param dt: time step
		 * @param T: final time
		 * @param Nx: number of grid points in x-direction
		 * @param Ny: number of grid points in y-direction
		 * @param a: reaction rate
		 * @param b: diffusion rate
		 * @param mu1: reaction term
		 * @param mu2: diffusion term
		 * @param eps: epsilon
		 */
		void setParameters(double dt, double T, int Nx, int Ny, double a, double b, double mu1, double mu2, double eps);


		/**
		* @brief Set the initial conditions for U and V
		* @details For y>ly/2, U = 1, 0 everywhere else
		* @details For x<lx/2, V = a/2, 0 everywhere else
		*/
		void setInitialConditions();

		/// @brief Solution loop
		void TimeIntegrate();

		/// @brief Write the solution to a file
		void writeToFile();

	private:
		/**
		 * @brief Solve f1
		 * @details Solve the reaction term f1 at a given node
		 * @param &u: pointer to the current U value
		 * @param &v: pointer to the current V value
		 */
		double solve_f1(double& u, double& v);

		/**
		 * @brief Solve f2
		 * @details Solve the diffusion term f2 at a given node
		 * @param &u: pointer to the current U value
		 * @param &v: pointer to the current V value
		 */
		double solve_f2(double& u, double& v);

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

		// U and V are double matrices of size Nx x Ny - represent the solution of the PDEs
		double *U, *V;

		// Derivative of U and V
		double *dU, *dV;
};


// Helper functions
void printFullMatrix(double* A, int Nx, int Ny);
void GNUPlot();

#endif