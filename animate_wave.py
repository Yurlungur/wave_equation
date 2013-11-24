#!/usr/bin/env python2

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-11-24 02:58:49 (jonah)>

# This is a companion program to my wave equation numerical solver. It
# generates movies from wave equation data.

# Imports
# ----------------------------------------------------------------------
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import animation
import sys,os
# ----------------------------------------------------------------------

# The number of types of field on the grid.
NUM_VAR_TYPES = 3 

def load_data(filename):
    """
    Takes a file name as a string and extracts the animation data from
    it. Returns a list of arrays:
    [times,energies,positions,frames]

    times and energies are one-dimensional arrays. They are the time
    and energy of each frame.

    frames is different. It is a three-dimensional array. The first
    index index the time of the frame. So frames[0] is the first
    frame. The second index is lattice point. So frames[0,0] is the
    left-hand side of the box and frames[0,-1] is the right-hand
    side. The third index indexes fields, ordered as: s,r,u.
    """
    with open(filename,'r') as f:
        lattice_spacing = eval(f.readline())
    raw_data = np.loadtxt(filename,skiprows=2)
    times = raw_data[...,0]
    energies = raw_data[...,-1]
    data = raw_data[...,1:-1]
    num_points = data.shape[-1]/NUM_VAR_TYPES
    num_frames = data.shape[0]
    frames = np.reshape(data,(num_frames,NUM_VAR_TYPES,num_points))
    positions = np.array([lattice_spacing * i for i in range(num_points)])
    return times,energies,positions,frames

def animate_data(times,positions,frames,filename):
    """
    Animates the data
    """
    fig = plt.figure()
    ax = plt.axes(xlim=(positions[0],positions[-1]),
                  ylim=(np.min(frames[0,-1],np.max(frames[0,-1]))))
    splot, = ax.plot([],[],'g-',label="s")
    rplot, = ax.plot([],[],'r-',label="r")
    uplot, = ax.plot([],[],'b-',label="u")
    time_template = 'time = %.1fs'
    time_text = ax.text(0.1*positions[-1],
                        0.9*np.max(frames[0,-1]),
                        '',transform = ax.transAxes)
    leg = ax.legend([splot,rplot,uplot]
                    ,["s","r","u"],loc=1)

    def init():
        splot.set_data([],[])
        rplot.set_data([],[])
        uplot.set_data([],[])
        time_text.set_text('')
        return splot,rplot,uplot,time_text

    def animate(i):
        splot.set_data(positions,frames[i,0])
        rplot.set_data(positions,frames[i,1])
        uplot.set_data(positions,frames[i,2])
        time_text.set_text(time_template%(times[i]))
        return splot,rplot,uplot,time_text

    anim = animation.FuncAnimation(fig,animate,frames=range(len(times)),
                                  init_func=init,blit=True)
    anim.save(filename+".mp4",fps=30)
    plt.show()

def plot_energy(times,energies):
    """
    Plots the energy as a function of time
    """
    plt.plot(times,energies)
    plt.xlabel('time')
    plt.ylabel(r'$r^2 + s^2$')
    plt.show()

def main(filename):
    times,energies,positions,frames = load_data(filename)
    plot_energy(times,energies)
    animate_data(times,positions,frames,filename)

for filename in sys.argv[1:]:
    main(filename)
    
