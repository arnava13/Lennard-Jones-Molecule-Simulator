# Diatomic Molecules Simulation

This repository contains a Python-based simulation of diatomic molecules' oscillations, developed for the Computer Modelling module at the University of Edinburgh's School of Physics and Astronomy. The simulation utilizes **Symplectic Euler** and **Velocity Verlet** integration algorithms to model vibrational behavior and evaluate the accuracy of the two methods.

## Introduction

The project simulates the oscillatory motion of diatomic molecules, specifically Nitrogen (N₂) and Oxygen (O₂), using numerical integration techniques. The report evaluates:
- Time interval parameters
- Energy conservation
- Wavenumber accuracy compared to experimental results

## Algorithms

### Symplectic Euler
- A first-order integration method that calculates position and velocity iteratively using Taylor expansions.
- Neglects higher-order terms, leading to noticeable energy conservation errors.

### Velocity Verlet
- A second-order integration method that calculates position and velocity using the time-averaged acceleration.
- Offers better accuracy and energy conservation.

## Physics Model

The simulation models the behavior of a diatomic molecule near its equilibrium separation:
- **Morse Potential**: Used to calculate forces between atoms.
- The oscillatory motion is approximated as simple harmonic:

$$ \ddot{r} = -\frac{D \alpha^2}{m}(r - r_e) $$

- Separation as a function of time:

$$ r(t) = r_0 \sin{\left(\sqrt{\frac{2 \alpha^2 D}{m}} t + \phi \right)} + r_e $$

- Total energy is used as a metric to assess accuracy:
  
$$ \text{Fractional Error} = \frac{E_{\text{max}} - E_{\text{min}}}{E_{\text{initial}}} $$

- Wavenumber is derived as:
  
$$ \tilde{\nu} = \frac{v}{c} $$

## Implementation

The simulation determines the optimal time interval by iterating over various values, targeting a maximum energy error of **5%**. The algorithm used to optimise the timestep is provided. Final interval:
- Chosen interval: **0.13 Å√(amu/eV)**
- Lower intervals were not used to reduce computational overhead.

## Input and Output Files

### Input Files
Input files are provided in **XYZ format (.dat)** and specify the initial conditions for spun and unspun molecules of Nitrogen (N₂) and Oxygen (O₂). These include:
- Atom type
- Initial positions
- Initial velocities
Similarly constructed XYZ input .dat files for other diatomic systems may be used.

### Output Files
The simulation writes results to an **XYZ format** file, including:
- Final positions of the molecules
- Final velocities

## Results

### Key Findings
- **Waveform**: Both algorithms produced sinusoidal oscillations, consistent with theoretical predictions.
- **Wavenumber**: Experimental vs simulated values show close agreement, though accuracy decreases for rotating molecules.
- **Energy Conservation**: Verlet significantly outperformed Euler in maintaining energy consistency.

### Results Table

| System Modelled | Algorithm | Wavenumber (cm⁻¹) | Error (%) | Energy Error (%) |
|-----------------|-----------|-------------------|-----------|------------------|
| Rotating O₂     | Euler     | 1383              | 12.4      | 3.1              |
| Rotating O₂     | Verlet    | 1389              | 12.1      | 0.1              |
| Non-Rotating O₂ | Euler     | 1530              | 0.3       | 3.9              |
| Non-Rotating O₂ | Verlet    | 1536              | 0.3       | 0.4              |
| Rotating N₂     | Euler     | 2203              | 6.6       | 3.7              |
| Rotating N₂     | Verlet    | 2216              | 6.1       | 0.4              |
| Non-Rotating N₂ | Euler     | 2309              | 2.1       | 5.2              |
| Non-Rotating N₂ | Verlet    | 2324              | 1.4       | 0.7              |

## Usage

1. **Clone the repository**:
 git clone https://github.com/yourusername/diatomic-molecules-simulation.git
 cd diatomic-molecules-simulation

2. Prepare input data: Ensure you have one of the provided input files in XYZ format. Inputs for the spun and unspun nitrogen and oxygen molecules are given.

3. Run the simulation: Use the following command to execute the program:

python exercise_3.py <input_file> <euler/verlet> <output_file>

<input_file>: Path to the input data file.
<euler or verlet>: Algorithm to use for the simulation (choose either euler or verlet).
<output_file>: Path to save the output data.

4. View the results

Check the <output_file> for final positions and velocities in XYZ format.
Visualise using a suitable .xyz file viewer.
