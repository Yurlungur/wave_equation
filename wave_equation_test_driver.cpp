// wave_equation_test_driver.cpp

// Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
// Time-stamp: <2013-11-23 23:40:03 (jonah)>

// This is my test driver for the wave equation simulation code.

// Includes
#include "wave_equation.hpp"
#include <iostream>
#include <fstream> // For file input and output
using namespace std;

int main() {
  dVector TEST_VECTOR;
  double lattice_spacing;
  ofstream output_file;
  string initial_data_filename = "initial_data.dat";

  cout << "----------------------------------------------\n"
       << "Test driver for the wave equation.\n"
       << "----------------------------------------------"
       << endl;

  cout << "\n\n----------------------------------------------\n"
       << "Testing the utility methods.\n"
       << "----------------------------------------------"
       << endl;
  cout << "Testing get_lattice_spacing." << endl;
  lattice_spacing = get_lattice_spacing(INTERVAL_LENGTH,NUM_POINTS);
  cout << "The lattice spacing is: " << lattice_spacing << endl;

  cout << "Testing the overloaded stream operator for dVectors."
       << endl;
  TEST_VECTOR.resize(4);
  for (int i = 0; i < 4; i++) {
    TEST_VECTOR[i] = i;
  }
  cout << "The test vector is: " << TEST_VECTOR << "." << endl;
  cout << "The test vector should be: [0, 1, 2, 3]." << endl;


  TEST_VECTOR.resize(3 * NUM_POINTS);
  for (int i = 0; i < 3 * NUM_POINTS; i++) {
    TEST_VECTOR[i] = i;
  }

  cout << "\n\n----------------------------------------------\n"
       << "Testing the array manipulation routines.\n"
       << "----------------------------------------------"
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
       << get_linear_index_u(50,TEST_VECTOR) << endl;
  cout << "The real index of the 99th element in r is: "
       << get_linear_index_r(99,TEST_VECTOR) << endl;
  cout << "The derivative of the zeroth element in s is: "
       << derivative(get_linear_index_s(0,TEST_VECTOR),
		     TEST_VECTOR,lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;
  cout << "The derivative of the 99th element of r is: "
       << derivative(get_linear_index_r(99,TEST_VECTOR),
		     TEST_VECTOR, lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;
  cout << "The derivative of the 50th element of u is: "
       << derivative(get_linear_index_u(50,TEST_VECTOR),
		     TEST_VECTOR,lattice_spacing)
       << ".\n"
       << "\tIt should be: " << float(1)/lattice_spacing << "." << endl;

  cout << "\n\n----------------------------------------------\n"
       << "Testing the initial data.\n"
       << "----------------------------------------------"
       << endl;
  cout << "Making initial standing wave with:\n"
       << "\tAmplitude = 1\n"
       << "\tWave number = 1\n"
       << "\tSpeed of sound = " << C2 << "\n"
       << "\tInterval length = " << INTERVAL_LENGTH << "\n"
       << "\tNum points = " << NUM_POINTS << endl;
  TEST_VECTOR = initial_standing_wave(1,1,C2,INTERVAL_LENGTH,NUM_POINTS);  
  cout << "The initial data is: " << endl;
  print_fields(TEST_VECTOR,lattice_spacing);
  cout << "\n" << endl;
  cout << "Printing to a the file: " << initial_data_filename << endl;
  // ofstreams want character arrays not strings, for some stupid
  // reason. The weird syntax feeds it the pointer to the first letter
  // in the string. Then it can figure it out from there.
  output_file.open(&initial_data_filename[0]);
  print_fields(TEST_VECTOR,lattice_spacing,output_file);
  output_file.close();
  cout << "Calculating the energy of the initial data.\n"
       << "\tEnergy is: " << energy(TEST_VECTOR,lattice_spacing) << "\n"
       << "\tEnergy should be: " << M_PI << endl;
  
  cout << "\n\n----------------------------------------------\n"
       << "This concludes the test.\n"
       << "----------------------------------------------"
       << endl;
}

