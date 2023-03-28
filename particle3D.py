"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Arnav Agarwal
Number: S2163065

"""
import numpy as np

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

    Attributes
    ----------
    label: name of the particle
    mass: mass of the particle
    position: position of the particle
    velocity: velocity of the particle

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
        """
        self.label = label
        self.mass = float(mass)
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.total_energy = 0 #To be updated by application.
        ...

    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        str
        """
        return self.label + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2])

    def kinetic_energy(self):
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        ke: float
            1/2 m v**2
        """
        v = np.linalg.norm(self.velocity)
        ke = 1/2*self.mass*v**2
        return ke
    
    def momentum(self):
        """
        Returns the momentum vector of a Particle3D instance
        
        Returns
        -------
        momentum: numpy.array
        """
        momentum = self.mass * self.velocity
        return momentum 
    
    def update_velocity(self, f, dt):
        """
        Updates the velocity of a Particle3D instance given a timestep and force.
        
        Parameters
        ----------
        
        f: numpy.array
            Force vector applied
        dt: float
            Timestep
        
        """
        self.velocity = self.velocity + f*dt/self.mass
        

    def update_position_1st(self, dt):
        """
        Updates the position of a Particle3D instance to the first order given a timestep.
        
        Parameters
        ----------
        dt: float
            Timestep
        
        """
        self.position = self.position + dt*self.velocity
        

    def update_position_2nd(self, f, dt):
        """
        Updates the position of a Particle3D instance to the second order given a timestep and a force.
        
        Parameters
        ----------
        f: numpy.array
            Force applied
            
        dt: float
            Timestep
        
        """
        self.position = self.position + dt*self.velocity + dt**2 * f * 1/(2*self.mass)

    @staticmethod
    def read_line(line):
        """
        Creates a Particle3D instance given a line of text.

        The input line should be in the format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

        Parameters
        ----------
        line: str
            Readable line in the above format

        Returns
        -------
        p: Particle3D
        """
        
        arguments = line.split()
        label = arguments[0]
        mass = float(arguments[1])
        position = [float(arguments[2]), float(arguments[3]), float(arguments[4])]
        velocity = [float(arguments[5]), float(arguments[6]), float(arguments[7])]
        p = Particle3D(label, mass, position, velocity)
        
        return p

    @staticmethod
    def total_kinetic_energy(particles):
        """
        Returns the total kinetic energy of a list of Particle3D instances
        
        Parameters
        ----------
        particles: list
                A list of Particle3D instances
                
        Returns
        ----------
        ke: float
            Total kinetic energy
        """
        
        ke = 0
        
        for particle in particles:
            ke += particle.kinetic_energy()
        
        return ke

    @staticmethod
    def com_velocity(particles):
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity
        """
        com_momentum = 0 
        com_vel = 0 
        total_mass = 0
        
        for particle in particles:
            com_momentum += particle.momentum()
            total_mass += particle.mass
            com_vel = com_momentum / total_mass
        return  com_vel
    



    
