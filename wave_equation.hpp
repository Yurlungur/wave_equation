// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-20 17:28:16 (jonah)>

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

// Access a variable at lattice point i of one type in the input
// dVector. 0 for s, 1 for r, and 2 for u.
double get_ith_element_of_type(int type, int i, const dVector& v);

// ----------------------------------------------------------------------
