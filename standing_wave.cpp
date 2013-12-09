// standing_wave.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-12-09 18:07:45 (jonah)>

// This is one piece of my wave equation numerical solver. This
// program is the main function for the simulation that looks for a
// standing wave solution.
// ----------------------------------------------------------------------


// Includes
#include "wave_equation.hpp"
#include <cassert> // for obvious reasons.
#include <cmath> // for math
#include <fstream> // For file input and output
using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::setw;


int main() {
  cout << "Let's simulate a standing wave!" << endl;
  
  // Values and objects we'll need
  std::string outfile_name = "standing_wave.dat"; // file name for data
  double interval_length = INTERVAL_LENGTH; // The length of the box
					    // the wave is confined to
  int num_points = NUM_POINTS; // The size of the grid
  int initial_data_algorithm = STANDING_WAVE; // initial data
  double boundary_conditions = OPEN; // Boundary conditions
  double lattice_spacing = get_lattice_spacing(interval_length,num_points,
					       boundary_conditions);
  double c2 = C2; // c^2, the speed of the wave squared.
  double max_t = 100*(2*M_PI); // The time we integrate to.
  double max_dt = 0.1*lattice_spacing/c2; // The max time step allowed
  double dt0 = max_dt/2.0; // The initial time step
  RKF45 integrator; // The runge kutta integrator object 

  // Let's do this!
  cout << "Working..." << endl;
  initialize_simulation(interval_length,num_points,c2,
			initial_data_algorithm,max_dt,dt0,
			boundary_conditions,integrator);
  integrator.integrate(max_t);
  cout << "Done!" << endl;

  // Find the energy change
  print_energy_statistics(integrator,lattice_spacing,boundary_conditions);

  // And let's print the output
  cout << "Printing...." << endl;
  print_animation_file(integrator,outfile_name,lattice_spacing,
		       boundary_conditions);
  cout << "Printing finished. Have a nice day." << endl;

  return 0;
}
