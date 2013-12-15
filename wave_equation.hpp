// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-12-11 17:34:34 (jonah)>

// This is the prototype for the wave equation program that uses the
// method of lines and finite differences to solve the wave equation.


// Background
// ----------------------------------------------------------------------
// See the included pdf or latex file.
// ----------------------------------------------------------------------


// Dependencies
// ----------------------------------------------------------------------
// This program depends on the RKF45 package I wrote, which is an
// implementation of a fourth-order Runge-Kutte Fehlberg integration
// scheme. If you don't have it, you can find it here:
// https://github.com/Yurlungur/runge_kutta
// ----------------------------------------------------------------------


// Usage
// ----------------------------------------------------------------------
// STUB
// ----------------------------------------------------------------------


// Include guard
#pragma once
// Includes
#include "rkf45.hpp" // The Runge-Kutta integrator
#include <vector> // for output and internal variables
#include <cmath> // for math
#include <iostream> // For input and output
#include <iomanip> // For input and output
#include <cstdlib>  // For random integers

// We use the dVector type a lot, so let's define a type for it
// to make things more readable.
typedef std::vector<double> dVector;

// ----------------------------------------------------------------------
// Global constants
// ----------------------------------------------------------------------

// Number of variable types. (i.e., s,r, and u. The integration variables.)
const int NUM_VAR_TYPES = 3;
// Debugging? 
const bool DEBUGGING = false;

// For now we will define the number of grid points and the interval
// length here. However, in a future modification, these will be moved
// to a configuration file.

// The spatial interval is [0,INTERVAL_LENGTH]
const double INTERVAL_LENGTH = 2*M_PI;
// The number of discrete grid points on the lattice
const int NUM_POINTS = 50;
// The Speed of the Wave
const double C2 = 1; // c^2

// Boundary data
const double LEFT_BOUNDARY_U = 0;  // left-hand boundary
const double RIGHT_BOUNDARY_U = 0; // right-hand boundary

// Some convenient names
const int STANDING_WAVE = 0; // for standing wave initial data
const int TRAVELLING_WAVE = 1; // For travelling wave initial data
const int S = 0; // For type s fields
const int R = 1; // For type r fields
const int U = 2; // for type u fields
// Is the boundary periodic?
const double PERIODIC = -1;
const double OPEN = 1;

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------

// Calculates the distance between points on the spatial grid.
double get_lattice_spacing(double interval_length, int num_points,
			   double boundary_conditions);

// Finds the difference between the largest element and the smallest
// element of a dVector. Useful for calculating energy differences.
double get_max_difference(const dVector& v);

// Find the max value of a dVector
double get_max(const dVector& v);

// Find hte minimum value of a dVector
double get_min(const dVector& v);

// Find the average value of a dVector
double get_average(const dVector& v);

// Find the variance of a dVector
double get_variance(const dVector& v);

// Get a random double between dmin and dmax. The distribution is
// uniform.
double uniform_distribution(double dmin,double dmax);

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Vector manipulation
// ----------------------------------------------------------------------

// The three variable vectors s,r, and u are all stored together in a
// large dVector. These methods allow easy access to them. We
// integrate u along with the other variables, although it doesn't
// contribute to the coupled system.

// Each variable is discretized in space to form a vector, indexed 0
// through NUM_POINTS-1. 

// The variables are all stored in a single extremely long
// dVector. The variables are stored in order: s,r,u. So the first
// NUM_POINTS elements are the elements of the s vector along the
// grid. And so fourth.

// Returns the number of points on the grid given a vector v. Used for
// internal calculations so the lattice size doesn't have to be passed
// around.
int get_num_points(const dVector& v);

// Access the real index of an element of type at lattice point i.
int get_linear_index_of_type(int type, int i, const dVector& v);

// Gets the real index of an element of s at lattice point i. 
int get_linear_index_s(int i, const dVector& v);

// Gets the real index of an element of r at lattice point i.
int get_linear_index_r(int i, const dVector& v);

// Gets the real index of an element of u at lattice point i.
int get_linear_index_u(int i, const dVector& v);

// Access a variable at lattice point i of one type in the input
// dVector. 0 for s, 1 for r, and 2 for u.
double get_ith_element_of_type(int type, int i, const dVector& v);

// Access the ith element of s
double get_ith_s(int i, const dVector& v);

// Access the ith element of r
double get_ith_r(int i, const dVector& v);

// Access the ith element of u
double get_ith_u(int i, const dVector& v);

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Calculus
// ----------------------------------------------------------------------

// Calculus operations on the grid are defined here. We define
// derivatives, inner products, and a notion of energy.

// Centered difference derivative. Takes the index i of the element
// and the dvector the element is in. Also needs the lattice
// spacing. This method has no error checking. Use at your own risk.
double centered_difference(int i, const dVector& v, double lattice_spacing);

// Centered difference derivative. Takes the index i of the element
// and the dvector the element is in. Also needs the lattice
// spacing. Is clever and uses modulo arithmetic to find the
// appropriate points for a periodic boundary.
double modular_centered_difference(int i, const dVector& v,
				   double lattice_spacing);

// Forward difference derivative. Takes the index i of the element and
// the dVector the element is in. Also needs the lattice spacing. This
// method has no error checking. Use at your own risk.
double forward_difference(int i, const dVector& v, double lattice_spacing);

// Backward difference derivative. Takes the index i of the element and
// the dVector the element is in. Also needs the lattice spacing. This
// method has no error checking. Use at your own risk.
double backward_difference(int i, const dVector& v, double lattice_spacing);

// Derivative operator that doesn't care about forward or backward
// differences. It automatically figures out whether or not the user
// is on the boundary and how to handle forward, centered, and
// backward differences. Takes the index i of the element, the vector
// the element is in, and the lattice spacing.
double derivative(int i, const dVector& v, double lattice_spacing,
		  double boundary_conditions);

// Derivative operator that acts on an element of r, s, or u. Otherwise
// like derivative. type 0 = s, type 1 = r, type 2 = u.
double derivative_for_type(int type, int i, const dVector& v,
			   double lattice_spacing,
			   double boundary_conditions);

// Derivative operator that acts on an element of s.
double derivative_for_s(int i, const dVector& v, double lattice_spacing,
			double boundary_conditions);

// Derivative operator that acts on an element of r.
double derivative_for_r(int i, const dVector& v, double lattice_spacing,
			double boundary_conditions);

// Derivative operator that acts on an element of u.
double derivative_for_u(int i, const dVector& v, double lattice_spacing,
			double boundary_conditions);

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Iteration Function
// ----------------------------------------------------------------------

// Here we define the iteration function of the form y' = f(t,y) which
// we will plug into the runge-kutta time stepper. f(t,y) now depends
// on spatial derivatives in y, which are just couplings between the
// elements of the y vector. The optional arguments vector contains
// the lattice spacing.
dVector f(double t, const dVector& y, const dVector& optional_args);

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Initial and boundary data
// ----------------------------------------------------------------------

// We need to define the initial conditions and the boundary
// conditions.  Since u(t,0) = alpha and u(t,L) = beta for all time, 
// s = (du/dt) is identically zero on the boundary.

// This in turn implies that (ds/dt) = 0 at the boundary for all
// time. We need to make sure that the initial data is consistent with
// this. So we must choose initial data such that
// (dr/dx) = 1/c^2 (ds/dt) is zero on the boundary.

// Makes initial data for a standing wave. Takes as input the
// amplitude of thew wave u0, the wave number k, the square speed of
// the wave c^2, the length of the interval L, and the number of
// points in the discretized grid, num_points.
// In the future we have many methods that generate initial data on
// the boundary. For now, we choose to make a standing wave with
// frequency and wave numbers of 1 and amplitude of 1/2.
// Function we want: u(x,t) = 2 u0 * cos(omega*t) * sin(k*x),
// where omega = k*c.
// Thus we want:
// u(0,x) = 2 * u_0 * sin(k * x)
// r(0,x) = 2 * u_0 cos(k * x)
// s(0,x) = 0 everywhere.
// amplitude = 2*u_0
dVector initial_standing_wave(double amplitude, double wave_number,
			      double c2, double interval_length,
			      int num_points);

// Makes initial data for a travelling wave. Takes as input the
// amplitude of the wave, the wave number k, the square speed of the
// wave c^2, the length of the interval L, and the number of points in
// the discretized grid, num_points.
// Function we want: u(x,t) = u0 * sin(kx - omega*t),
// where omega = k*c.
// Thus we want:
// u(0,x) = u0*sin(kx)
// r(0,x) = k*u0*cos(kx)
// s(0,x) = -omega*u0*cos(kx)
dVector initial_travellng_wave(double amplitude, double wave_number,
			       double c2, double interval_length,
			       int num_points);

// Makes initial data for the system. A helper
// method. initial_data_type determines the initial data function to
// call.
dVector get_initial_data(double amplitude, double wave_number,
			 double c2, double interval_length,
			 int num_points,
			 int initial_data_type);
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Energy
// ----------------------------------------------------------------------

// We need a way to measure the energy of the system. Ideally, this
// does not grow. The energy is
// E = integral(r^2 + s^2) over the grid.
// To make this work, we need a notion of inner product. See the
// background for more information.

// Inner product between two vector types on the grid, type_1 and
// type_2. These are types of vector field as defined above. Requires
// the lattice spacing of the grid.
double inner_product(int type_1, int type_2, const dVector& grid,
		     double lattice_spacing,double boundary_conditions);

// Calculates the energy of the grid. Requires the lattice spacing.
double energy(const dVector& grid, double lattice_spacing,
	      double boundary_conditions);

// Get the energy at each time step of a finished integrator
// object. Useful for bug checking.
dVector get_all_energies(const RKF45& integrator, double lattice_spacing,
			 double boundary_conditions);

// Get the the total change in energy over the course of a simulation
// from the integrator
double get_energy_difference(const RKF45& integrator,
			     double lattice_spacing,
			     double boundary_conditions);

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Input and output
// ----------------------------------------------------------------------

// Methods for generating data.

// Prints a snapshot of the fields on the input grid.
void print_fields(const dVector& grid, double lattice_spacing,
		  std::ostream& out = std::cout);

// Overload the stream input operator to print an input dVector. If
// you print the fields on a grid, this will look pretty bad.
std::ostream& operator<< (std::ostream& out, const dVector& in);

// Prints the data data for an integrated animation. The format is one
// line per frame, with each line containing:
// time s_field r_field u_field energy
// where the fields are all space-separated lists. Energy and time are
// the energy of the system and time of the frame
// respectively. lattice_spacing is the distance between points on the
// grid. The first line contains the lattice spacing.
void print_animation_data(const RKF45& system, std::ostream& out,
			  double lattice_spacing,
			  double boundary_conditions);

// Prints the animation file. Defined as above.
void print_animation_file(const RKF45& system, const std::string& filename,
			  double lattice_spacing,
			  double boundary_conditions);

// Prints some statistics about the energy of a finished system as a
// function of time.
void print_energy_statistics(const RKF45& system,
			     std::ostream& out,
			     double lattice_spacing,
			     double boundary_conditions);

// Prints some statistics about the energy of a finished system as a
// function of time.
void print_energy_statistics(const RKF45& system, double lattice_spacing,
			     double boundary_conditions);
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Simulation
// ----------------------------------------------------------------------

// We generate initial data and feed it into the Runge-Kutta
// algorithm. Input the type of initiald data you want (currently the
// only option is standing wave, which is type 0.

// Interval length is the length of the interval. num_points is the
// number of grid points. c2 is the square speed of the
// wave. initial_data_algorithm is the algorithm used to generate the
// initial data. max_dt is the time step size. Should not be chosen
// dynamically if you wish to make a movie. dt0 is the initial step
// size. Play with it to get a good value. integrator is the
// Runge-Kutta integrator used. It is assumed to be fresh with nothing
// set and no time steps made. After the algorithm is done, the
// integrator will be returned.
void initialize_simulation(double interval_length, int num_points,
			   double c2,
			   int initial_data_algorithm,
			   double max_dt, double dt0,
			   double boundary_conditions,
			   RKF45& integrator);
// Default initial data algorithm is STANDING_WAVE
void initialize_simulation(double interval_length, int num_points,
			   double c2,
			   double max_dt, double dt0,
			   RKF45& integrator);

// After all is said and do

// ----------------------------------------------------------------------
