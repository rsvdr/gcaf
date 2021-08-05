# gcaf

This guiding center code can simulate test ions in a maxwellian plasma composed of electrons and a single species of ions subject to an external model magnetic field corresponding to a tokamak.

## Software structure

* **main.f90** Contains the time loop of the simulation.
* **shared.f90** Contains shared constants and variables among subroutines.
* **initial.f90** Contains the initialization of variables.
* **step.f90** Contains the equations of motion that evolve the particle variables.
* **collisions.f90** Contains the Boozer and Kuo-Petravic collision operators.
* **read.f90** Contains subroutines that read the user defined parameters.
* **record.f90** Contains subroutines that write the .plt files.
* **functions.f90** Contains functions used among subroutines.
* **field.f90** Contains the analytic magnetic field and the field perturbation.

## Installation

You can install GCAF simply with the command
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
number_of_particles      4
time                     800      in on-axis toroidal transits
step_delta               2000     toroidal transits/N
step_fix                 off      on/off
max_error                1e-8     maximum energy error in each step
trajectory_points        1000     number of points in traj.plt
dispersion_points        400      number of points in dips.plt
distribution_points      4        number of distribution snapshots
poincare_lower_bound     0.4      in r/rmin
poincare_upper_bound     0.6      in r/rmin

# FIELD
field_magnitude          2        in 10kG
major_radius             150      in cm
minor_radius             50       in cm
perturbation_amplitude   5e-4
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
density                  2e14     in cm^-3
density_profile          flat     flat/para

# TEST PARTICLE
particle_type            ion      ion/ele
mass                     1        in proton units
charge                   1        in electron charge units
pitch_angle              0        pitch angle in grades
energy                   1        in keV
energy_distribution      mon      mon/max
radius                   -1       in cm
theta                    0        in grades
zeta                     0        in grades

# OUTPUT FILES
traj.plt    tim rad tht zet par per eng pch rcy zcy
poin.plt    rad tht rcy zcy
dist.plt    rad eng pch
disp.plt    tim dsr dse dsz skr ske skz flr fle flz
coef.plt    den flx dfr dfe dfp sdr sde sdp
```

One can run a simulation simply as ```gcaf -i input_directory/data.in -o output_directory/```.

Some command-line flags that can be used are:
```
-sim  Type of simulation traj/poin/ense
-col  Collisions on/off
-prt  Perturbation on/off
-npr  Number of particles
-pch  Initial pitch angle (in degrees)
-eng  Initial energy in keV
-amp  Field perturbation amplitude
-pot  Electrostatic potential in keV
-den  Plasma density in cm^-3
```
