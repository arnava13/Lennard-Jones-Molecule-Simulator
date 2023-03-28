#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Symplectic Euler and Velocity Verlet time integration of
a particle moving in a morse potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to an output file.

The potential is V(r) = D((1-exp(-a(r-r_e)))^2-1), where
D, a and r_e are obtained from the input file and passed to the functions that
calculate force and potential energy.

Author: Arnav Agarwal
Number: S2163065.

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from scipy import signal
from particle3D import Particle3D

def morse_force(D, a, r_e, p1, p2):
    """
    Return the magnitude of the force on particle 1 due to particle 2
    in a morse potential.

    The force is given by
        F(r) = 2aD((1-exp(-a(r-r_e)))*exp(-a(r-r_e))).

    Parameters
    ----------
    p1: Particle3D
        Particle3D instance for particle 1.
    p2: Particle3D 
        Particle3D instance for particle 2.
    D: float 
        Depth of potential well.
    a: float 
        Curvature of potential well.
    r_e: float 
        Equilibrium separation.

    Returns
    -------
    force: numpy.ndarray
    """
    r = np.linalg.norm(p1.position - p2.position)
    r_12 = p2.position - p1.position 
    rhat_12 = r_12 / np.linalg.norm(r_12)
    force = 2*a*D*((1-np.exp(-a*(r-r_e)))*np.exp(-a*(r-r_e)))
    return force * rhat_12


def morse_potential(D, a, r_e, p1, p2):
    """
    Returns the potential energy of particles in a morse potential,
    given by
    V(r) = D((1-exp(-a(r-r_e)))^2-1).

    Parameters
    ----------
    p1: Particle3D
        Particle3D instance for particle 1.
    p2: Particle3D 
        Particle3D instance for particle 2.
    D: float 
        Depth of potential well.
    a: float 
        Curvature of potential well.
    r_e: float 
        Equilibrium separation.
    
    Returns
    -------
    potential: float
    """
    r = np.linalg.norm(p1.position - p2.position)
    potential = D*((1-np.exp(-a*(r-r_e)))**2-1)
    return potential

def period_wavenumber_frequency(separations, timeunit, dt):
    """
    Finds the average period of oscillations and converts it from microscopic units to seconds.
    Uses this.
    Uses this average to estimate the wavenumber of emitted photons from the system, and converts this value from
    inverse-metres to inverse-centimetres.
    Also uses average period to estimate the natural frequency of the particles' oscillations, in Hertz (inverse-seconds).
    

    Parameters
    ----------
    separations : 1D Array-like
        Input array describing seperation of two particles at the ith time instant.

    Returns
    -------
    period : float
        Average period of oscillations in seconds.
    wavenumber: float
        Estimated wavenumber of emitted photons in cm^-1.
    frequency: float
        Estimated natural frequency of oscillations in Hertz.
        
        

    """
    #Find positions of peaks using Scipy.signal.
    peaks = signal.find_peaks(separations)
    peaktimes = peaks[0]
    
    #Sum over all measured periods (difference in indices of subsequent peaks),
    #and divide by number of measured periods to find average period.
    periods_measured = 0
    periods_summed = 0
    for i in range(1, len(peaktimes)):
        period = peaktimes[i] - peaktimes[i-1]
        periods_summed += period
        periods_measured += 1
    period_weirdunits = periods_summed/periods_measured
    
    #Convert measured period from Angstrom / sqrt(AMU/ Electronvolt) to seconds.
    period_seconds = period_weirdunits * timeunit * dt
    
    #Calculate natural frequency of the oscillations.
    frequency = 1/period_seconds
    
    #Find wavenumber from the period (=1/Tc=vc) where v is the natural frequency, c is the speed of light, and convert to cm^-1.
    c = 299792458
    wavenumber_permetre = frequency/c
    wavenumber_percentimetre = wavenumber_permetre/100
    
    return period_seconds, wavenumber_percentimetre, frequency


def main():
    # Read inputs from command line
    # The variable sys.argv contains whatever you typed on the command line
    # or after %run on the ipython console to launch the code.  We can use
    # it to get user inputs.
    # Here we expect three things:
    #    the name of this file
    #    the name of the input data file
    #    euler or verlet
    #    the name of the output file the user wants to write to
    # So we start by checking that all three are specified and quit if not,
    # after giving the user a helpful error message.
    if len(sys.argv) != 4 :
        print("You did not provide exactly three arguments.")
        print("In spyder, run like this instead:")
        print(f"    %run {sys.argv[0]} <desired input data file> <euler or verlet> <desired output file>")
        sys.exit(1)
    else:
        infile_name = sys.argv[1]
        mode = sys.argv[2]
        outfile_name = sys.argv[3]
        
    # Open the input file for reading ("r"), and read its lines.
    infile = open(infile_name, "r")
    infile_lines = infile.readlines()

    # Open the output file for writing ("w")
    outfile = open(outfile_name, "w")
    
    #Set starting time at zero.
    time = 0.0

    # Read simulation parameters from input file.
    sim_parameters = infile_lines[1].split()
    dt = float(sim_parameters[0])
    numstep = int(sim_parameters[1])
    # numstep = 100/dt
    #I used the commented line above in my report in order to keep the simulation
    #time constant.
    
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

    # Get magnitude of initial force
    force = morse_force(D, a, r_e, p1, p2)
    
    #Calculate initial force vector on each particle
    f1 = force
    f2 = -force #By Newton's third law.

    # Write out starting time, position, and energy values
    # to the output file.
    energy = p1.kinetic_energy() + p2.kinetic_energy() + morse_potential(D, a, r_e, p1, p2)
    outfile.write(f"{time}    {p1.position}    {p2.position}    {energy}\n")

    # Initialise numpy arrays that we will plot later, to record
    # the trajectories of the particles.
    times = np.zeros(numstep)
    separations = np.zeros(numstep)
    energies = np.zeros(numstep)
    
    #Append initial energy and separations to respective arrays:
    energies[0] = energy
    separations[0] = np.linalg.norm(p1.position - p2.position)

   
    for i in range(1, numstep):

        # Update the positions and velocities.
        # This will depend on whether we are doing an Euler or verlet integration
        if mode == "euler":
            # Update particle position to 1st order
            p1.update_position_1st(dt)
            p2.update_position_1st(dt)
            
            # Get the force at new positions.
            f1 = morse_force(D, a, r_e, p1, p2)
            f2 = -f1 #By newton's third law.
            

            # Update particle velocity 
            p1.update_velocity(f1, dt)
            p2.update_velocity(f2, dt)

        elif mode == "verlet":
        
            # Update particle position to 2nd order using previous force
            p1.update_position_2nd(f1, dt)
            p2.update_position_2nd(f2, dt)
            
            # Get the force at the new positions.
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

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + p2.kinetic_energy() + morse_potential(D, a, r_e, p1, p2)
        outfile.write(f"{time}    {p1.position}    {p2.position}    {energy}\n")

        # Store the things we want to plot later in our arrays
        times[i] = time
        separations[i] = np.linalg.norm(p1.position - p2.position)
        energies[i] = energy
        
        

    # Now the simulation has finished we can close our output file
    outfile.close()

    # Plot particle trajectory to screen. There are no units
    # here because it is an abstract simulation, but you should
    # include units in your plot labels!
    modelabel = mode[0].upper() + mode[1:]
    pyplot.figure()
    pyplot.title(modelabel + ': Unspun $O_2$ Particle Separation vs Time')
    pyplot.xlabel('Time $(\AA\sqrt{amu/eV})$')
    pyplot.ylabel('Particle Separation $(\AA)$')
    ax = pyplot.gca()
    ax.ticklabel_format(useOffset=False)
    pyplot.plot(times, separations)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.figure()
    pyplot.title(modelabel + ': Total Energy vs Time')
    pyplot.xlabel('Time $(\AA\sqrt{amu/eV})$')
    pyplot.ylabel('Energy $(eV)$')
    ax = pyplot.gca()
    ax.ticklabel_format(useOffset=False)
    pyplot.plot(times, energies)
    pyplot.show()
    
    #Print the estimated frequency of the oscillations and the estimated wavenumber of photons.
    timeunit = 1.018050571 * 10**-14
    oscillation_outcomes = period_wavenumber_frequency(separations, timeunit, dt)
    print(f"\nThe frequency of oscillations of the two-particle system is approximately {oscillation_outcomes[2]} Hz.\n")
    print(f"The wavenumber of emitted photons is approximately {oscillation_outcomes[1]} cm^-1.\n")
    
    #Find the deviation in maximum and minimum energy of oscillator.
    delta_E = np.abs(np.max(energies) - np.min(energies))
    #Find and print the energy innacuracy as delta_E divided by the initial energy.
    init_energy = np.abs(energies[0])
    energy_innacuracy = delta_E/init_energy
    print(f"The ratio of the variation in calculated maximum and minimum total energy to the initial total energy is {energy_innacuracy}.\n")
    
    
    


# This strange but standard python idiom means that the main function will
# only be run if we run this file, not if we just import it from another
# python file. It is good practice to include it whenever your code can be
# run as a program.
if __name__ == "__main__":
    main()

