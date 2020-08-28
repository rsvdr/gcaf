# gcafpp

This guiding center code can simulate test ions in a maxwellian plasma composed of electrons and a single species of ions subject to an external model magnetic field corresponding to a tokamak.

For more detailed information about the code please visit the [wiki page](https://github.com/Saavestro/gcaf/wiki/GCAF-wiki).

## Software structure

Some contents in the source files are the following:

* **main.cpp** Contains the particle and time loop of the simulation.
* **shared.h** Contains the shared constants and variables among subroutines.
* **define.h** Contains definitions.
* **initial.cpp** Initializes simulation variables.
* **step.cpp** Contains the equations of motion that evolve particle variables.
* **collisions.cpp** Contains the Boozer and Kuo-Petravic collision operators.
* **read.cpp** Contains subroutines that read user defined parameters.
* **record.cpp** Contains subroutines that record data and write the .plt files.
* **functions.cpp** Contains functions used among subroutines.
* **field.cpp** Contains the analytic magnetic field and the field perturbation.

## Prerequisites

You need gcc and Python3

## Installation

You can install gcafpp with the command
```
make && make install
```
This will produce an executable file that can be run from the command-line if you include ```~/.local/bin``` in your ```$PATH```.

## Usage

A sample input file is

```
# SWITCHES
simulation               traj     traj/ense/poin
device                   tkmk     tkmk/stell
collisions               off      on/off
perturbation             off      on/off

# SIMULATION
number_of_particles      200
time                     80       in on-axis toroidal transits
step_delta               2000     toroidal transits/N
step_fix                 off      on/off
max_error                1e-8     maximum energy error in each step
trajectory_points        4000     number of points in traj.plt
moment_points            4000     number of points in dips.plt
distribution_points      4        number of distribution snapshots
poincare_lower_bound     0.1      in r/rmin
poincare_upper_bound     0.9      in r/rmin

# FIELD
field_magnitude          2        in 10kG
major_radius             150      in cm
minor_radius             50       in cm
perturbation_amplitude   1e-4
perturbation_frequency   0        in Hz
perturbation_phase       0        in grades
m_mode                   2
n_mode                   1
electrostatic_potential  0        in keV
potential_profile        cons     pres/cons

# PLASMA
ions_mass                1        in proton units
ions_charge              1        in electron charge units
ions_temperature         1        in keV
electrons_temperature    1        in keV
temperature_profile      flat     flat/gaus
density                  2e15     in cm^-3
density_profile          flat     flat/para

# PARTICLE
particle_type            ion      ion/ele
mass                     1        in proton units
charge                   1        in electron charge units
pitch_angle              0        pitch angle in grades
energy                   1.0      in keV
energy_distribution      mon      mon/max
radius                   25       in cm
theta                    45       in grades
zeta                     0        in grades
```

One can run a simulation as ```gcafpp -i input_directory/data.in -o output_directory/```.

Some command-line flags that can be used are:
```
-sim   Type of simulation traj/poin/ense
-npr   Number of particles
-prt   Field Perturbation on/off
-amp   Field perturbation amplitude
-col   Collisions on/off
-den   Plasma density in cm^-3
-pch   Initial pitch angle (in degrees)
-eng   Initial energy in keV
-pot    Electrostatic potential in keV
```

## Postprocesing with Python

The output .plt data files produced in a simulation can be postprocessed using Python. Various postprocessing programs can be found in the folder pys. These can be run as `./file.py -i input_directory -o output_directory`.

* **traj.py** Plot particle trajectories and the time evolution of particle variables obtained in a trajectory simulation.
* **dist.py** Plot a snapshot of the particle, energy and pitch angle distribution function obtained in an ensemble simulation.
* **mome.py** Plot of the time evolution of the statistical moments obtained in an ensemble simulation.
* **poin.py** Plot a Poincaré section.
* **surf.py** Postprocess a given particle Poincaré surface. The output data is needed for a flux calculation.
