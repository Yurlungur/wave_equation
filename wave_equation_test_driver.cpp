// wave_equation_test_driver.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-20 17:31:36 (jonah)>

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

  cout << "Testing the array manipulation routines." << endl;
  cout << "The fiftieth element in s is: "
       << get_ith_element_of_type(0,50,TEST_VECTOR) << endl;
  cout << "The fiftieth element in r is: "
       << get_ith_element_of_type(1,50,TEST_VECTOR) << endl;
  cout << "The fiftieth element in u is: "
       << get_ith_element_of_type(2,50,TEST_VECTOR) << endl;
  
  cout << "This concludes the test." << endl;
}
