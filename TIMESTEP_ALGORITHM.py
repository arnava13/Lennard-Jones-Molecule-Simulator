#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Timestep optimisation algorithm
for Symplectic Euler and Velocity Verlet 
time integration of a particle moving in a morse potential.

Produces a plot of the fractional error in total energy versus timestep in 
Å*sqrt(amu/eV).

Author: Arnav Agarwal
Number: S2163065.

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from scipy import signal
from exercise_3 import morse_force, morse_potential
from particle3D import Particle3D

def main():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect three things:
    #    the name of this file
    #    euler or verlet
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) != 3 :
        print("You did not provide exactly three arguments.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <desired input data file> <euler or verlet>")
        sys.exit(1)
    else:
        infile_name = sys.argv[1]
        mode = sys.argv[2]
        
    # Open the input file for reading ("r"), and read its lines.
    infile = open(infile_name, "r")
    infile_lines = infile.readlines()
    
    # Determine range of timesteps to be plotted.
    dts = np.linspace(0.1, 0.2, 10)
    
    
    # Read morse potential parameters from input file.
    morse_parameters = infile_lines[3].split()
    D = float(morse_parameters[0])
    a = float(morse_parameters[2])
    r_e = float(morse_parameters[1])

    # Read Particle3D parameters for two particles from input file: label, mass, position, velocity.
    # Initialise two Particle3D instances based on these parameters.
    p1_parameters = infile_lines[4]
    p2_parameters = infile_lines[5]
    p1 = Particle3D.read_line(p1_parameters)
    p2 = Particle3D.read_line(p2_parameters)

    #Initialise arrays to be filled with the energy error values.
    energy_errors = np.zeros(len(dts))
    
    #Initialise index for looping through dts.
    j=-1

    # Loop through various values of dt. 
    for dt in dts: 
        #Increase dt loop index by 1.
        j += 1
        
        # Get initial forces.
        f1 = morse_force(D, a, r_e, p1, p2)
        f2 = -f1 #By Newton's third law.
        
        #Determine number of timesteps needed for total of 100[T] simulation time.
        numstep = int(100/dt)
        
        # Initialise numpy arrays energy of the particle at each instant.
        energies = np.zeros(numstep)
        
        # Calculate Initial Energy
        energy_init = p1.kinetic_energy() + p2.kinetic_energy() + morse_potential(D, a, r_e, p1, p2)
        
        #Append initial energy to energies arrays:    
        energies[0] = energy_init
        
        #Start the time integration loop.
        for i in range(1, numstep):
            # Update the positions and velocities.
            # This will depend on whether we are doing an Euler or verlet integration
            if mode == "euler":
                # Update particle positions
                p1.update_position_1st(dt)
                p2.update_position_1st(dt)
                
                # Calculate force at new positions
                f1 = morse_force(D, a, r_e, p1, p2)
                f2 = -f1 #By newton's third law.
                
    
                # Update particle velocity 
                p1.update_velocity(f1, dt)
                p2.update_velocity(f2, dt)
    
            elif mode == "verlet":
            
                # Update particle position using previous force
                p1.update_position_2nd(f1, dt)
                p2.update_position_2nd(f2, dt)
                
                # Get the force for the new positions.
                f1_new = morse_force(D, a, r_e, p1, p2)
                f2_new = -f1_new #By newton's third law.
                
                # Update particle velocity by averaging
                # current and new forces
                p1.update_velocity(0.5*(f1+f1_new), dt)
                p2.update_velocity(0.5*(f2+f2_new), dt)
                
                # Re-define force value for the next iteration
                f1 = f1_new
                f2 = f2_new
            else:
                raise ValueError(f"Unknown mode {mode} - should be euler or verlet")
    
            # Store the things we need for our error calculations.
            energy = p1.kinetic_energy() + p2.kinetic_energy() + morse_potential(D, a, r_e, p1, p2)
            energies[i] = energy
        
        #Find energy error for this dt value.
        delta_E = np.max(energies) - np.min(energies)
        energy_error = np.abs(delta_E/energy_init)
    
        #Make this error the jth value in the energy errors array.
        energy_errors[j] = energy_error
    
    # Plot patricle energy error vs timestep to screen.
    modelabel = mode[0].upper() + mode[1:]
    pyplot.figure()
    pyplot.title(modelabel + ': Fractional Error in Total Energy vs Timestep')
    pyplot.xlabel('Timestep $(\AA\sqrt{amu/eV})$')
    pyplot.ylabel('Fractional error in Total Energy')
    ax = pyplot.gca()
    ax.ticklabel_format(useOffset=False)
    pyplot.plot(dts, energy_errors)
    pyplot.show()


# This strange but standard python idiom means that the main function will
# only be run if we run this file, not if we just import it from another
# python file. It is good practice to include it whenever your code can be
# run as a program.
if __name__ == "__main__":
    main()

