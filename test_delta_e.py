#!/usr/bin/env python2

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-11-27 11:45:20 (jonah)>

# This is a simple script that takes the energy data from a .dat file
# generated by the wave equation simulator and finds the total change
# in energy over the course of the simulation. Call it on a single (or
# many files) with
# python2 test_delta_e.py standing_wave.dat

# Imports
# ----------------------------------------------------------------------
import numpy as np
import sys,os
import animate_wave # Another python program with data loading commands
from decimal import Decimal # a way to get nice formatting for float error
# ----------------------------------------------------------------------

def get_energy_difference(filename):
    """
    Gets the energy for a simulation
    """
    times,energies,positions,frames = animate_wave.load_data(filename)
    delta = np.max(energies) - np.min(energies)
    return delta

if __name__ == "__main__":
    for filename in sys.argv[1:]:
        print "{}".format(Decimal(get_energy_difference(filename)))