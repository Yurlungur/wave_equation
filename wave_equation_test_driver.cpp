// wave_equation_test_driver.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-21 10:43:06 (jonah)>

// This is my test driver for the wave equation simulation code.

// Includes
#include "wave_equation.hpp"
#include <iostream>
using namespace std;

int main() {
  dVector TEST_VECTOR;
  double lattice_spacing;

  cout << "Testing the uitility methods." << endl;
  cout << "Testing get_lattice_spacing." << endl;
  lattice_spacing = get_lattice_spacing(INTERVAL_LENGTH,NUM_POINTS);
  cout << "The lattice spacing is: " << lattice_spacing << endl;

  for (int i = 0; i < 3 * NUM_POINTS; i++) {
    TEST_VECTOR.push_back(i);
  }

  cout << "----------------------------------------------\n"
       << "Testing the array manipulation routines.\n"
       << "----------------------------------------------\n"
       << endl;
  cout << "The number of points in the array is: "
       << get_num_points(TEST_VECTOR) << endl;
  cout << "The fiftieth element in s is: "
       << get_ith_element_of_type(0,50,TEST_VECTOR) << endl;
  cout << "The fiftieth element in r is: "
       << get_ith_element_of_type(1,50,TEST_VECTOR) << endl;
  cout << "The fiftieth element in u i1s: "
       << get_ith_element_of_type(2,50,TEST_VECTOR) << endl;
  cout << "The real index of the fiftieth element in u is: "
       << get_real_index_u(50,TEST_VECTOR) << endl;
  cout << "The real index of the 99th element in r is: "
       << get_real_index_r(99,TEST_VECTOR) << endl;
  cout << "The derivative of the zeroth element in s is: "
       << derivative(get_real_index_s(0,TEST_VECTOR),
		     TEST_VECTOR,lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;
  cout << "The derivative of the 99th element of r is: "
       << derivative(get_real_index_r(99,TEST_VECTOR),
		     TEST_VECTOR, lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;
  cout << "The derivative of the 50th element of u is: "
       << derivative(get_real_index_u(50,TEST_VECTOR),
		     TEST_VECTOR,lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;
  cout << "This concludes the test." << endl;
}
