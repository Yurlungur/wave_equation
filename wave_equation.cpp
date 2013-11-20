// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-20 16:58:21 (jonah)>

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
#include <iostream> // For input and output


// ----------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------
// Calculates the distance between points on the spatial grid. 
double get_lattice_spacing(double interval_length, int num_points) {
  return double(interval_length)/(num_points - 1);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Vector manipulation
// ----------------------------------------------------------------------
// Access a variable at lattice point i of one type in the input
// dVector. 0 for s, 1 for r, and 2 for u.
double get_ith_element_of_type(int type, int i, const dVector& v) {
  assert ( 0 <= type && NUM_VAR_TYPES > type
	   && "The type must be one of the allowed variable types.");
  assert ( v.size() % NUM_VAR_TYPES == 0
	   && "The vector must be for the variables of the system.");
  // A local "num points" variable. Used for simplicity.
  int num_points = v.size() / 3; // Guaranteed to be an even integer
				 // by above assert statement.
  assert (0 <= 1 && i < num_points
	  && "The ith elemennt must be on the grid.");
  
  // The true index of the element we want in the array.
  int element_index = type * num_points + i;

  return v[element_index];
}
