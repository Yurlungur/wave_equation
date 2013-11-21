// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-21 11:18:43 (jonah)>

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
#include <vector> // for output and internal variables
#include <cmath> // for math

// We use the dVector type a lot, so let's define a type for it
// to make things more readable.
typedef std::vector<double> dVector;


// ----------------------------------------------------------------------
// Global constants
// ----------------------------------------------------------------------

// Number of variable types. (i.e., s,r, and u. The integration variables.)
const int NUM_VAR_TYPES = 3;

// For now we will define the number of grid points and the interval
// length here. However, in a future modification, these will be moved
// to a configuration file.

// The spatial interval is [0,INTERVAL_LENGTH]
const double INTERVAL_LENGTH = 2*M_PI;
// The number of discrete grid points on the lattice
const int NUM_POINTS = 100;
// The Speed of the Wave
const double C2 = 1; // c^2
// Debugging? 
const bool DEBUGGING = false;

// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------

// Calculates the distance between points on the spatial grid.
double get_lattice_spacing(double interval_length, int num_points);

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
int get_real_index_of_type(int type, int i, const dVector& v);

// Gets the real index of an element of s at lattice point i. 
int get_real_index_s(int i, const dVector& v);

// Gets the real index of an element of r at lattice point i.
int get_real_index_r(int i, const dVector& v);

// Gets the real index of an element of u at lattice point i.
int get_real_index_u(int i, const dVector& v);

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
double derivative(int i, const dVector& v, double lattice_spacing);

// Derivative operator that acts on an element of r, s, or u. Otherwise
// like derivative. type 0 = s, type 1 = r, type 2 = u.
double derivative_for_type(int type, int i, const dVector& v,
			   double lattice_spacing);

// Derivative operator that acts on an element of s.
double derivative_for_s(int i, const dVector& v, double lattice_spacing);

// Derivative operator that acts on an element of r.
double derivative_for_r(int i, const dVector& v, double lattice_spacing);

// Derivative operator that acts on an element of u.
double derivative_for_u(int i, const dVector& v, double lattice_spacing);

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
