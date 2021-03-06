/*      The Reaction-Diffusion Problem (HPC-fied)!
        CID: 01707131       */

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <omp.h>

#include "ReactionDiffusion.h"

namespace po = boost::program_options;


int main(int argc, const char *argv[]){

    // Boost program options
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help", "produce help message")
        ("threads", po::value<int>()->default_value(16), "Run the solver with specified number of threads.")
        ("dt", po::value<double>()->default_value(0.001), "Time-step to use.")
        ("T", po::value<double>()->default_value(100.0), "Total integration time.")
        ("Nx", po::value<int>()->default_value(101), "Number of grid points in x.")
        ("Ny", po::value<int>()->default_value(101), "Number of grid points in y.")
        ("a", po::value<double>()->default_value(0.75), "Value of parameter a.")
        ("b", po::value<double>()->default_value(0.06), "Value of parameter b.")
        ("mu1", po::value<double>()->default_value(5.0), "Value of parameter mu1.")
        ("mu2", po::value<double>()->default_value(0.0), "Value of parameter mu2.")
        ("eps", po::value<double>()->default_value(50.0), "Value of parameter epsilon.");

    // Default input for execution will look like:
    // ./main --threads 16 --dt 0.001 --T 10.0 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --mu1 5.0 --mu2 0.0 --eps 50.0

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << "\n";
        return 1;
    }

    // Loading parameters into memory
    const double dt     = vm["dt"].as<double>();
    const double T      = vm["T"].as<double>();
    const int Nx        = vm["Nx"].as<int>();
    const int Ny        = vm["Ny"].as<int>();
    const double a      = vm["a"].as<double>();
    const double b      = vm["b"].as<double>();
    const double mu1    = vm["mu1"].as<double>();
    const double mu2    = vm["mu2"].as<double>();
    const double eps    = vm["eps"].as<double>();

    // Set the number of threads
    omp_set_num_threads(vm["threads"].as<int>());


    std::cout << "\nWelcome to the Reaction-Diffusion Solver! Brace yourself for epic speeds...\n" << std::endl;


    // Initializing the solver class
    ReactionDiffusion myAwesomeSolver;

    // Setting the parameters
    myAwesomeSolver.setParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);

    // Setting the initial conditions
    myAwesomeSolver.setInitialConditions();

    // Solving the PDEs!
    myAwesomeSolver.TimeIntegrate();

    // Writing the solution to output.txt
    myAwesomeSolver.writeToFile();


    std::cout << "\nBye!" << std::endl;

    return 0;

}