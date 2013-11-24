// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-24 00:45:18 (jonah)>

// This is the implementation of the wave equation program that uses
// the method of lines and finite differences to solve the wave
// equation.


// For background, see the included pdf or latex file.
// For usage, see the readme or the header file.


// Dependencies
// ----------------------------------------------------------------------
// This program depends on the RKF45 package I wrote, which is an
// implementation of a fourth-order Runge-Kutte Fehlberg integration
// scheme. If you don't have it, you can find it here:
// https://github.com/Yurlungur/runge_kutta
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


// ----------------------------------------------------------------------
double get_lattice_spacing(double interval_length, int num_points) {
  return double(interval_length)/(num_points - 1);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
int get_num_points(const dVector& v) {
  assert ( v.size() % NUM_VAR_TYPES == 0
	   && "The vector must be for the variables of the system.");
  return v.size() / 3;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
int get_linear_index_of_type(int type, int i, const dVector& v) {
  assert ( 0 <= type && NUM_VAR_TYPES > type
	   && "The type must be one of the allowed variable types.");
  // A local "num points" variable. Used for simplicity.
  int num_points = get_num_points(v);
  assert (0 <= i && i < num_points
	  && "The ith elemennt must be on the grid.");
  
  // The true index of the element we want in the array.
  int element_index = type * num_points + i;

  return element_index;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
int get_linear_index_s(int i, const dVector& v) {
  return get_linear_index_of_type(S,i,v);
}

int get_linear_index_r(int i, const dVector& v) {
  return get_linear_index_of_type(R,i,v);
}

int get_linear_index_u(int i, const dVector& v) {
  return get_linear_index_of_type(U,i,v);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
double get_ith_element_of_type(int type, int i, const dVector& v) {
  return v[get_linear_index_of_type(type,i,v)];
}

double get_ith_s(int i, const dVector& v) {
  return get_ith_element_of_type(S,i,v);
}

double get_ith_r(int i, const dVector& v) {
  return get_ith_element_of_type(R,i,v);
}

double get_ith_u(int i, const dVector& v) {
  return get_ith_element_of_type(U,i,v);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
double centered_difference(int i, const dVector& v, double lattice_spacing) {
  return (v[i+1] - v[i-1])/(2*lattice_spacing);
}

double forward_difference(int i, const dVector& v, double lattice_spacing) {
  return (v[i+1] - v[i])/lattice_spacing;
}

double backward_difference(int i, const dVector& v, double lattice_spacing) {
  return (v[i] - v[i-1])/lattice_spacing;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
double derivative(int i, const dVector& v, double lattice_spacing) {
  // A local "num points" variable. Used for simplicity.
  int num_points = get_num_points(v);
  assert (0 <= i && i < (int)v.size()
	  && "The ith elemennt must be in the vector");

  if ( i % num_points == 0 ) { // We're at the lower boundary of the grid.
    if (DEBUGGING) {
      cout << "\tUsing forward difference." << endl;
    }
    return forward_difference(i,v,lattice_spacing);
  }
  else if (i % num_points == num_points - 1) {
    // We're at the upper boundary of the grid.
    if (DEBUGGING) {
      cout << "\tUsing backward difference." << endl;
    }
    return backward_difference(i,v,lattice_spacing);
  }
  else { // We're in the middle of the grid.
    if (DEBUGGING) {
      cout << "\tUsing centered difference." << endl;
    }
    return centered_difference(i,v,lattice_spacing);
  }
}
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
double derivative_for_type(int type, int i, const dVector& v,
			   double lattice_spacing) {
  return derivative(get_linear_index_of_type(type,i,v),v,lattice_spacing);
}

double derivative_for_s(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(S,i,v,lattice_spacing);
}

double derivative_for_r(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(R,i,v,lattice_spacing);
}

double derivative_for_u(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(U,i,v,lattice_spacing);
}
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Iterative Function
// ----------------------------------------------------------------------
dVector f(double t, const dVector& y, const dVector& optional_args) {
  // t is not used. But it is required by the integrator.

  // The optional arguments vector contains the lattice spacing and c^2.
  double h = optional_args[0];  // h = lattice spacing
  double c2 = optional_args[1]; // c2 = c^2.

  // For convenience, get the number of points in the lattice
  int num_points = get_num_points(y);

  // The output dVector. Starts empty.
  dVector output(y.size(),0);

  // The order of elements in v is s,r,u. So we have the following system:
  // (d/dt) [s,r,u] = [c^2 (dr/dx), (ds/dx), s]

  // Fill the output array
  for (int i = 0; i < num_points; i++) {
    // Fill the (du/dt) terms with s.
    output[get_linear_index_u(i,output)] = get_ith_s(i,y);
    // Fill the (dr/dt) vectors with (ds/dx).
    output[get_linear_index_r(i,output)] = derivative_for_s(i,y,h);
    // Fill the (ds/dt) vectors with c^2(dr/dx).
    output[get_linear_index_s(i,output)] = c2 * derivative_for_r(i,y,h);
  }
  // (du/dt) at the boundary is forced to be zero. Impose this.
  output[get_linear_index_s(0,output)] = 0;
  output[get_linear_index_s(num_points-1,output)] = 0;

  // That's it. Return the array!
  return output;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Boundary data
// ----------------------------------------------------------------------
dVector initial_standing_wave(double amplitude, double wave_number,
			      double c2, double interval_length,
			      int num_points) {
  assert ( c2 > 0 && "The square speed of the wave must be positive." );

  // The first thing we need to do is calculate the lattice spacing
  double lattice_spacing = get_lattice_spacing(interval_length,num_points);

  // We also need the output:
  // A vector with the appropriate number of elements for the three
  // fields on the lattice.
  dVector initial_data(NUM_VAR_TYPES * num_points,0);

  // A temporary variables for storing u, r, and s at a given lattice
  // point.
  double u,r,s;

  // s is initially zero everywhere.
  s = 0;

  // Now we need to fill the array with the initial data we care about
  for (int i = 0; i < num_points; i++) {
    u = amplitude * sin(wave_number * lattice_spacing * i);
    r = amplitude * cos(wave_number * lattice_spacing * i);
    initial_data[get_linear_index_u(i,initial_data)] = u;
    initial_data[get_linear_index_r(i,initial_data)] = r;
    initial_data[get_linear_index_s(i,initial_data)] = s;
  }
 
  // Return the generated standing wave.
  return initial_data;  
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void print_fields(const dVector& grid, double lattice_spacing,
		  std::ostream& out) {
  // Header line
  int precision = 10;
  int width = precision+5;
  out << std::setiosflags(std::ios::left)
      << std::setfill(' ')
      << setw(width) << "#x "
      << setw(width) << " s"
      << setw(width) << "  r"
      << setw(width) << "   u" << endl;

  // Formatting goes here if necessary
  out << std::setprecision(precision)
      << std::setiosflags(std::ios::left);
  //      << std::setfill('0');

  // Get the number of grid points
  int num_points = get_num_points(grid);

  // Output!
  for (int i = 0; i < num_points; i++) {
    out << setw(width) << i*lattice_spacing << " "
	<< setw(width) << grid[get_linear_index_s(i,grid)] << " "
	<< setw(width) << grid[get_linear_index_r(i,grid)] << " "
	<< setw(width) << grid[get_linear_index_u(i,grid)] << " "
	<< endl;
  }
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void print_animation_data(const RKF45& system, std::ostream& out,
			  double lattice_spacing) {
  // Print header line
  out << lattice_spacing << endl;
  out << "# s_field r_field u_field energy time" << endl;

  // For each time step, we will need a y vector. We'll re-use this
  // one.
  dVector y;
  // We'll need a time too. Use this one.
  double t;
  // And an energy. Use this one.
  double e;

  // Iterate through every time step of the system and output our data.
  for (int n = 0; n < system.steps(); n++) {
    t = system.get_t(n);
    y = system.get_y(n);
    e = energy(y,lattice_spacing);
    out << t << " ";
    for (dVector::const_iterator it = y.begin(); it != y.end(); ++it) {
      out << (*it) << " ";
    }
    out << e << endl;
  }
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void print_animation_file(const RKF45& system, const std::string& filename,
			  double lattice_spacing) {
  ofstream outfile;
  outfile.open(&filename[0]);
  print_animation_data(system,outfile,lattice_spacing);
  outfile.close();
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
std::ostream& operator <<(std::ostream& out, const dVector& in) {
  out << "[";
  for (unsigned int i = 0; i < in.size()-1; i++) {
    cout << in[i] << ", ";
  }
  cout << in[in.size()-1] << "]";
  return out;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Energy
// ----------------------------------------------------------------------
double inner_product(int type_1, int type_2, const dVector& grid,
		     double lattice_spacing) {
  // For convenience, the number of points on the grid
  int num_points = get_num_points(grid);
  // The output number
  double output = 0;
  // For convenience, the product of just two elements
  double product;

  for (int i = 0; i < num_points; i++) {
    product = get_ith_element_of_type(type_1,i,grid)
      * get_ith_element_of_type(type_2,i,grid);
    if ( i == 0 || i == num_points - 1 ) { // metric has 1/2s here
      output += 0.5 * lattice_spacing * product;
  }
    else { // metric is the identity here.
      output += lattice_spacing * product;
    }
  }
  return output;
}

double energy(const dVector& grid, double lattice_spacing) {
  return inner_product(0,0,grid,lattice_spacing)
    + inner_product(1,1,grid,lattice_spacing);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Simulation
// ----------------------------------------------------------------------
void initialize_simulation(double interval_length, int num_points,
			   double c2,
			   int initial_data_algorithm,
			   double max_dt, double dt0,
			   RKF45& integrator) {

  // First generate initial data.
  assert (initial_data_algorithm == STANDING_WAVE
	  && "There are no other algorithms available at this time.");
  // We assume the following values for amplitude and wave number
  double amplitude = 1;
  double wave_number = 1;
  dVector y0 = initial_standing_wave(amplitude,wave_number,c2,
				     interval_length,num_points);

  // We need to calculate the lattice spacing.
  double lattice_spacing = get_lattice_spacing(interval_length,num_points);

  // Our iterator function takes optional arguments
  dVector optional_args(2);
  optional_args[0] = lattice_spacing;
  optional_args[1] = c2; // The square speed of the wave.

  // Let's set the other important values for the integrator.
  integrator.set_f(f); // Iterator function
  integrator.set_y0(y0); // initial data
  integrator.set_t0(0); /// initial time
  integrator.set_optional_args(optional_args); // fed into the iterator function
  integrator.set_max_dt(max_dt); // Max time step
  integrator.set_dt0(dt0); // initial step size 
}

void initialize_simulation(double interval_length, int num_points,
			   double c2,
			   double max_dt, double dt0,
			   RKF45& integrator) {
  initialize_simulation(interval_length,num_points,c2,
			STANDING_WAVE,max_dt,dt0,integrator);
}
// ----------------------------------------------------------------------
