#include "ForwardEuler.h"

using namespace std;

// Constructor
ForwardEuler::ForwardEuler(){
	// Open output file
	output.open( "output.txt", ios::out | ios::trunc );    
}



void Solver::xy_conv(){
	// Convert degrees of freedom v_1 and v_2 to 2D cartesian coordinates
	x[0] = l[0] * sin(v[0]);
	y[0] = -l[0] * cos(v[0]);
	x[1] = x[0] + l[1] * sin(v[1]);
	y[1] = y[0] - l[1] * cos(v[1]);
}

void Solver::write(){

	// Write calculated coordinates to 'output.txt'
	output.precision(6);
	output << setw(15) << (iters * h)
		   << setw(15) << x[0] 
		   << setw(15) << y[0] 
		   << setw(15) << x[1] 
		   << setw(15) << y[1] << endl;
}