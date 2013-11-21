// wave_equation.hpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-21 11:18:43 (jonah)>

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
using std::cout;
using std::endl;


// ----------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
double get_lattice_spacing(double interval_length, int num_points) {
  return double(interval_length)/(num_points - 1);
}
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Vector manipulation
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
int get_num_points(const dVector& v) {
  assert ( v.size() % NUM_VAR_TYPES == 0
	   && "The vector must be for the variables of the system.");
  return v.size() / 3;
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
int get_real_index_of_type(int type, int i, const dVector& v) {
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
int get_real_index_s(int i, const dVector& v) {
  return get_real_index_of_type(0,i,v);
}

int get_real_index_r(int i, const dVector& v) {
  return get_real_index_of_type(1,i,v);
}

int get_real_index_u(int i, const dVector& v) {
  return get_real_index_of_type(2,i,v);
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
double get_ith_element_of_type(int type, int i, const dVector& v) {
  return v[get_real_index_of_type(type,i,v)];
}

double get_ith_s(int i, const dVector& v) {
  return get_ith_element_of_type(0,i,v);
}

double get_ith_r(int i, const dVector& v) {
  return get_ith_element_of_type(1,i,v);
}

double get_ith_u(int i, const dVector& v) {
  return get_ith_element_of_type(2,i,v);
}
// ----------------------------------------------------------------------



// ----------------------------------------------------------------------
// Calculus
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
  return derivative(get_real_index_of_type(type,i,v),v,lattice_spacing);
}

double derivative_for_s(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(0,i,v,lattice_spacing);
}

double derivative_for_r(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(1,i,v,lattice_spacing);
}

double derivative_for_u(int i, const dVector& v, double lattice_spacing) {
  return derivative_for_type(2,i,v,lattice_spacing);
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
    output[get_real_index_u(i,output)] = get_ith_s(i,y);
    // Fill the (dr/dt) vectors with (ds/dx).
    output[get_real_index_r(i,output)] = derivative_for_s(i,y,h);
    // Fill the (ds/dt) vectors with c^2(dr/dx).
    output[get_real_index_s(i,output)] = c2 * derivative_for_r(i,y,h);
  }

  // That's it. Return the array!
  return output;
}
// ----------------------------------------------------------------------
