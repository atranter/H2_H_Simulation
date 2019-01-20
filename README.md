# Molecular Kinematic Simulator

Author: William Simon

## Purpose
Through the use of [OpenFermion](https://github.com/quantumlib/OpenFermion), 
[OpenFermion-Psi4](https://github.com/quantumlib/OpenFermion-Psi4),
[Psi4](http://psicode.org/), and Python, the files in this repository are used to simulate the 
kinematics of each atom in a molecule, specified by the user on input. The data generated from 
the `simulation.py` file can then be used as input to the `plot.py` to display the evolution
of the molecule. These animations can also be saved to mp4 format using the `ffmpeg` package.

With the increasing research surrounding Quantum Computing and its applications, this program
will be used for research intended to demonstrate the possible applications of Quantum Computing to
computational chemistry. 

Through the command line, the user is able to specify various molecular paramaters that increase the
modularity of this program such as:
  * Each Atom in the Simulation 
  * Initial Cartesian Coordinates of Each Atom
  * Initial Velocity of Each Atom in Each Axis
  * Multiplicity 
  * Charge 
  
Additionally, global variables in the `simulation.py` file can be changed to increase modularity.
These variables will be described in more detail below along with descriptions of the theory, workflow, 
and usage of this program.

## Theory
To estimate the forces acting on every atom during an iteration, each atom is displaced by a small amount
(in both positive and negative directions) in every coordinate axis. Next, for each configuration, OpenFermion
is used to calculate the ground state energy. The energy differences between the displaced configuration
and the original configuration are then used to estimate the average force acting on each atom in all three
coordinate axes.

To then project the velocity and location of each atom in every coordinate axis, two methods are used. The less
accurate method, called the Euler-Cromer method, will work for most molecules and can be calculated quicker, 
however the timestep for each iteration. The second method, called the Runge-Kutta 4th order method, can project
out to larger timesteps and has increased accuracy, but requires a significantly longer computation time. 

After each iteration, the current geometry and the current velocity of every atom is stored for the next iteration.

## Workflow
In order to realize the full intended usage of this program, the client should first use the 
`simulation.py` script to generate the location data of each atom in the input molecule. Then,
the `plot.py` script can be used to display the data generated in the form of animations depicting
each atom as either a point in space or a curve showing the trajectory of the molecule.

## Usage
**simulation.py**

`python simulation.py [output datafile] [multiplicity] [charge]`

Followed by:

`[atomic symbol] [z coordinate] [y coordinate] [x coordinate] [z velocity] [y velocity] [x velocity]`

for each atom intended to be in the simulation.

For example: `python simulation.py example.txt 1 0 H 0 0 0 0 0 0 H 0 0 0.96382 0 0 0` 
simulates the H2 molecule at a bond length of +30% from equilibrium and outputs the data to example.txt

**plot.py**

`python plot.py [datafile] [y/n]`

Where [y/n] (either yes or no) is an optional input and determines whether or not to save the 
created animation in mp4 format to disk.

For example: `python plot.py sim_data/H2_30%_1000.txt n` shows the evolution of the H2 molecule
with a bond length of +30% from equilibrium over 1000 iterations. To show the animation as either
points or lines, uncomment the associated line in the `plot.py` file under the main function. 

## TODO:
  * requirements.txt
  

NOTE: I do not take any ownership for the `PlotMotionPeter.py` file. 
