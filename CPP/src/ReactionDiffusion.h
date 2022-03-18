#ifndef REACTION_DIFFUSION_H
#define REACTION_DIFFUSION_H

#include <vector>
#include "ForwardEuler.h"

class ReactionDiffusion : public ForwardEuler{
	private:
		// These attributes are specific to the RK4 method hence made private
		std::vector<std::vector<double>> K;
		
	public:
		// Constructor
		ReactionDiffusion();
		void solve();
};


#endif