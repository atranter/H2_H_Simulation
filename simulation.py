''' 
Purpose: This script uses the OpenFermionWrapper (dependent on OpenFerimon and
Psi4) to simulate the kinematics of molecules given through command line inputs.
This simulation works by using small displacements in the location of each atom
in the z, y, and x axes and calculating the ground state energy of the
associated hamiltonian (produced by OpenFermion) to determine the force acting
on each atom in the molecule. With this determined force, either the
Euler-Cromer or Runge-Kutta 4th order methods are used to predict the velocity
and location of each atom at every point in the time of the simualtion. The
location of each atom is printed in a single line to be used in subsequent 
programs. The client can change the number of iterations and the length of
each timestep by setting the appropriate global variables below (ITERATIONS and
dt respectively)

INPUT FORMAT (all on one line):
			   python simulation.py [datafilename] [multiplicity] [charge]/
for each atom: [atomic symbol] [z coord] [y coord] [x coord]/
			   [velocity z] [velocity y] [velocity x]

EXAMPLE FORMAT:
python simulation.py example.txt 1 0 H 0 0 0 0 0 0 H 0 0 0.7414 0 0 0
^above simulates the H2 molecule at equilibrium bond-length with no initial
velocity for either atom in any axis

TODO:
	print hamil for RK4
	change RK4 to calculate all accelerations at the same time
	change Initial velocities and geometries from global variables to local
	change velocity data structure to resemble geometry
	make it easier to edit the geometry
	can we make anything multithreaded? - hard with psi4 calculations

Biggest Time Concern:
	calling get_ground_state() is by far the bottle-neck of this program

'''

from OpenFermionWrapper import OpenFermionWrapper
import sys
import math

import timeit

''' --- Set Constants ---'''
dt = 5*10**-8 #s from fs
dx = 0.001 #angstrom
dy = 0.001 #angstrom -> how much a molecule is displaced to find nearby energies
dz = 0.001 #angstrom

ITERATIONS = 300

BASIS = "sto-3g" #default basis... note -> a larger basis set may dramatically
				# increase runtime for solving for ground state energy

MASSES = list()

# first multiple is for MeV/c^2 to eV/c^2, second multiple is for amu to MeV
# data from https://www.lenntech.com/periodic/mass/atomic-mass.htm
MASS_DICT = {"H":  1000000*931.5*1.0079,
			 "He": 1000000*931.5*4.002602,
			 "Li": 1000000*931.5*6.941,
			 "Be": 1000000*931.5*9.0122,
			 "B":  1000000*931.5*10.811,
			 "C":  1000000*931.5*12.0107,
			 "N":  1000000*931.5*14.0067,
			 "O":  1000000*931.5*15.9994,
			 "F":  1000000*931.5*18.9984,
			 "Ne": 1000000*931.5*20.1797,
			 "Na": 1000000*931.5*22.9897,
			 "Mg": 1000000*931.5*24.305,
			 "Al": 1000000*931.5*26.9815,
			 "Si": 1000000*931.5*28.0855,
			 "P":  1000000*931.5*30.9738,
			 "S":  1000000*931.5*32.065,
			 "Cl": 1000000*931.5*35.453,
			 "K":  1000000*931.5*39.0983,
			 "Ar": 1000000*931.5*39.948,
			 "Ca": 1000000*931.5*40.078,
			 "Sc": 1000000*931.5*44.9559,
			 "Ti": 1000000*931.5*47.867,
			 "V":  1000000*931.5*50.9415,
			 "Cr": 1000000*931.5*51.9961,
			 "Mn": 1000000*931.5*54.938,
			 "Fe": 1000000*931.5*55.845,
			 "Ni": 1000000*931.5*58.6934,
			 "Co": 1000000*931.5*58.9332,
			 "Cu": 1000000*931.5*63.546,
			 "Zn": 1000000*931.5*65.39,
			 "Ga": 1000000*931.5*69.723,
			 "Ge": 1000000*931.5*72.64,
			 "As": 1000000*931.5*74.9216,
			 "Se": 1000000*931.5*78.96,
			 "Br": 1000000*931.5*79.904,
			 "Kr": 1000000*931.5*83.8}

''' --- From Command Line/Input --- '''
DATAFILE = "error.txt"
INITIAL_VELOCITIES = list()
INITIAL_GEOMETRY = list()
N = 0 #number of atoms in simulation
MULTIPLICITY = 0
CHARGE = 0

def parse_inputs(input):
	''' 
	ARGS: 
		input: command line parameters

	RETURNS:
		void

	This function parses the inputs and sets the global variables 
	appropriately. It also checks that the user has input the correct
	number of parameters. If an incorrect number of parameters are input,
	then the program exits. The number of atoms is determined by the length
	of the input and the program requires 7 inputs per atom (atomic symbol,
	z-coord, y-coord, x-coord, z-velocity, y-velocity, x-velocity. It then
	also obtains the multiplicity and charge of the molecule via the user.

	''' 

	global DATAFILE, MULTIPLICITY, CHARGE, N, dt
	global INITIAL_VELOCITIES, INITIAL_GEOMETRY, MASSES, MASS_DICT

	# each atom produces an additional 7 inputs and there should be 4 other
	# inputs, meaning there needs to be at least 11 and should be 7n+4
	if(((len(input)-4)%7 != 0) or (len(input) < 11)):
		sys.exit("Incorrect number of input parameters. See documention.\n\n")

	# set global variables
	DATAFILE = input[1]
	MULTIPLICITY = int(input[2])
	CHARGE = int(input[3])
	N = ((len(input)-4)/7)

	for i in range(N):
		# index of the current atom in the command line arguements
		atom_i = (i*7)+4

		# mass for each atom is precalculated and stored in MASS_DICT
		MASSES.append(MASS_DICT[str(input[atom_i])])

		# atoms are stored in the OpenFermion geometry syntax:
		# ([atomic symbol], ([z], [y], [x]))
		atom = tuple((str(input[atom_i]), (float(input[atom_i+1]), 
				   float(input[atom_i+2]), float(input[atom_i+3]))))
		INITIAL_GEOMETRY.append(atom)

		# velocities are stored as one array
		INITIAL_VELOCITIES.append(float(input[atom_i+4]))
		INITIAL_VELOCITIES.append(float(input[atom_i+5]))
		INITIAL_VELOCITIES.append(float(input[atom_i+6]))


def write_hamiltonians_to_file(hamiltonian):
	''' 
	ARGS: 
		input: qubit hamiltonian

	RETURNS:
		void

	This function casts the input qubit hamiltonian (typically from OpenFerimon)
	to a string and then prints it to a file. This is done at periodic intervals
	to show the evoltion of the hamiltonian over time.

	''' 
	filename = "hamil_" + DATAFILE
	file = open(filename, "a")
	file.write(str(hamiltonian))
	file.write("\n")
	file.write("\n")
	file.close()


def write_data(geometry):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))

	RETURNS:
		void

	This function prints the geometry of the current molecule to a data file 
	that shows the evolution of the atomic coordinates over time. The structure
	merely prints the z, y, and x coordinates (in order) for each atom 
	(sequentially) on the same line. It is expected that this function is called
	at every iteration so that the data can then be used in subsequent programs.

	''' 
	data = open(DATAFILE, "a")
	for i in range(N):
		if(i == N-1):
			# If the last atom, don't print a space at the end
			data.write("{} {} {}".format(geometry[i][1][0], geometry[i][1][1],
									     geometry[i][1][2]))
		else:
			data.write("{} {} {} ".format(geometry[i][1][0], geometry[i][1][1],
									      geometry[i][1][2]))
	data.write("\n")
	data.close()


def getGroundState(geometry):
    ''' 
    ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
	RETURNS:
		ground state energy: float in eV
		qubit_hamiltonian: generated hamiltonain from OpenFermion
    
    This function obtains the ground state energy of the given molecule
    	and the qubit hamiltonian using the OpenFermionWrapper class, the
    	OpenFermion program, and the Psi4 program. It then converts from units 
    	of Hartree and returns the energy of the molecule in eV. '''
    molecule = OpenFermionWrapper()
    molecule.load_molecule(geometry=geometry,         basis=BASIS, 
    					   multiplicity=MULTIPLICITY, charge=CHARGE,
    					   forceCalculation=True)
    molecule.set_ground_state_energy()
    molecule.perform_transform("BK")
    # Convert from Hartree to eV
    # return 27.2114*molecule.molecule.hf_energy
    return 27.2114*molecule.ground_state_energy, molecule.qubit_hamiltonian_bk


def displace_geometry(geometry, dz, dy, dx, a):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		dz: displacement in the z axis (float)
		dy: displacement in the y axis (float)
		dx: displacement in the x axis (float)
		a: index of the atom in the geometry list to displace

	RETURNS:
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))

	This function is designed to displace a single atom in the geometry given
	as input and return a geometry object (in the same format) with the intended
	atom displaced.
	''' 
	displaced_geometry = list()
	for i in range(N):
		if(i != a): # not the atom to be updated, so just append the current pos
			displaced_geometry.append(geometry[i]) 
		else: # atom to be slightly displaced
			displaced_geometry.append(
				tuple((str(geometry[i][0]), # atomic symbol
			    (float(geometry[i][1][0])+dz, 
			     float(geometry[i][1][1])+dy,     # update each coordinate
			     float(geometry[i][1][2])+dx))))
	return displaced_geometry
	
	
def get_time_step(vi, d, a):
	''' 
	ARGS: 
		vi: initial veloicity (float)
		d: distance to cover (float)
		a: acceleration (float)

	RETURNS:
		time: float

	This function calculates the time it would take for an atom with the current
	velocity and acceleration to traverse the distance provided as input. 
	''' 
	vf = math.sqrt((vi**2) + 2*a*d)
	time = (2*d)/(vi+vf)
	return time


def set_time_step(velocities, forces):
	''' 
	ARGS: 
		velocities: list of velocities of each atom in 3 axes (array of floats)
		forces: list of forces of each atom in 3 axes (array of floats)
	RETURNS:
		void

	This function attempts to set the time step for the next iteration by
	determining the time it would take for an atom to traverse a set distance.
	This method would be used for a variational timestep simulation.'''
	global dt
	# get the index of the atom with the largest acting force 
	# Note: this should maybe be velocities not forces
	index = forces.index(max(forces, key=abs))
	mass = MASSES[index/3]
	# calculate the acceleration on the given atom
	a = (forces[index]/mass)*10**10
	# set the global variable dt with the calculated time
	dt = get_time_step(velocities[index], 0.01, a)


def update_positions_EC(geometry, velocities):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		velocities: list of velocities of each atom in 3 axes (array of floats)
	RETURNS:
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))

	This function uses the Euler-Cromer method to update the positions of each
	atom after a time of dt using the velocities stored in the velocities array.
	''' 
	newGeometry = list()
	for i in range(N): # for each atom in our molecule, update the position
		z_coord = float(geometry[i][1][0]) + velocities[(i*3)]*dt
		y_coord = float(geometry[i][1][1]) + velocities[(i*3)+1]*dt
		x_coord = float(geometry[i][1][2]) + velocities[(i*3)+2]*dt

		newGeometry.append(tuple((str(geometry[i][0]), # atomic symbol
					             (z_coord, y_coord, x_coord)))) # position
	return newGeometry


def findForces_EC(geometry, hamil):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		hamil: a boolean indicating whether to save the hamiltonian to disk or
			ignore it
	RETURNS:
		average_forces: a list of floats describing the force acting on each
			atom in each of the 3 axes

	This function uses the Euler-Cromer method to update the positions of each
	atom after a time of dt using the velocities stored in the velocities array.
	''' 
	# create empty local lists to store local variables 
	geometries = list()
	energies = list()
	forces = list()

	# list of forces to return
	average_forces = list()

	geometries.append(geometry)

	# for each atom in the molecule, displace the atom in positive and negative
	# x and y axis and store the resulting geometry
	for i in range(N):
		# right now, only doing x and y directions, not z
		_r_geometry = displace_geometry(geometry, 0, 0, dx, i)
		geometries.append(_r_geometry)
		_l_geometry = displace_geometry(geometry, 0, 0, -dx, i)
		geometries.append(_l_geometry)
		_u_geometry = displace_geometry(geometry, 0, dy, 0, i)
		geometries.append(_u_geometry)
		_d_geometry = displace_geometry(geometry, 0, -dy, 0, i)
		geometries.append(_d_geometry)

	# for each of the stored geometries, calculate and store the energy
	for g in geometries:
		# Energy given in eV
		energy, hamiltonian = getGroundState(g)
		energies.append(energy)
		# write hamiltonian to disk if requested
		if(hamil):
			write_hamiltonians_to_file(hamiltonian)

	# for each of the stored energies, calculate the force acting on the atom
	# in all directions (positive and negative on each axis) by calculating the
	# difference in energy from current geometry to the intial
	for e in energies:
		# Store as eV/m
		forces.append(-1*(e-energies[0])/(dx*10**-10))

	# calculate and store the average force between the positive and negative
	# directions
	for i in range(N):
		# index is i*3+1?
		average_x = (forces[(4*i)+1]-forces[(4*i)+2])/2
		average_y = (forces[(4*i)+3]-forces[(4*i)+4])/2
		average_forces.append(0) # not worrying about forces in z yet
		average_forces.append(average_y)
		average_forces.append(average_x)

	return average_forces


def update_velocities_EC(velocities, forces):
	''' 
	ARGS: 
		velocities: list of velocities of each atom in 3 axes (array of floats)
		forces: list of forces of each atom in 3 axes (array of floats)
	RETURNS:
		velocities: list of velocities of each atom in 3 axes (array of floats)

	This function calculates the velocity of each atom in all 3 coordinate axes
	and returns the updated list. The velocities are calculated by evolving each
	atom (with respecitve mass) by the associated force (given as input). 
	''' 
	# if using variational timestep, uncomment line below
	# set_time_step(velocities, forces)

	# for each atom, calculate the new velocity in all 3 axes after the set time
	for i in range(0, N):
		velocity_z = velocities[(i*3)]
		velocity_y = velocities[(i*3)+1]
		velocity_x = velocities[(i*3)+2]
		fz = forces[(i*3)]
		fy = forces[(i*3)+1]
		fx = forces[(i*3)+2]
		velocity_z += (fz/MASSES[i])*dt*10**10
		velocity_y += (fy/MASSES[i])*dt*10**10
		velocity_x += (fx/MASSES[i])*dt*10**10
		velocities[(i*3)] = velocity_z
		velocities[(i*3)+1] = velocity_y
		velocities[(i*3)+2] = velocity_x
	return velocities


def calc_single_force_RK4(geometry, atom_i, axis):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		atom_i: int representing the index of the current atom
		axis: int representing the axis to act on (0 - z) (1 - y) (2 - x)
	RETURNS:
		float: average force acting on specified atom in specified axis

	This function calculates the force acting on a single atom in the desired
	axis. Similar to the find_forces_EC() function, this calculation is done
	by displacing the atom in the positive and negative directions by a fixed
	distance, calculating the associated energy of the hamiltonian, and
	averaging the difference between the positive and negative displacements
	with the energy of the current configuration.
	''' 
	geometries = list()
	energies   = list()
	forces     = list()

	#  Used to get energy of current location
	geometries.append(geometry)

	#  Get geometry of the molecule with one atom displaced along one axis
	if(axis == 0):
		geometries.append(displace_geometry(geometry,  dz,   0,   0, atom_i))
		geometries.append(displace_geometry(geometry, -dz,   0,   0, atom_i))
	elif(axis == 1):
		geometries.append(displace_geometry(geometry,   0,  dy,   0, atom_i))
		geometries.append(displace_geometry(geometry,   0, -dy,   0, atom_i))
	elif(axis == 2):
		geometries.append(displace_geometry(geometry,   0,   0,  dx, atom_i))
		geometries.append(displace_geometry(geometry,   0,   0, -dx, atom_i))

	for g in geometries:
		# Energy given in eV
		energy, hamiltonian = getGroundState(g)
		energies.append(energy)

	for e in energies:
		# Store as eV/m
		forces.append(-1*(e-energies[0])/(dx*10**-10))

	return (forces[1]-forces[2])/2


def calc_accel_RK4(vel, geometry, atom_i, axis, timestep, mass):
	''' 
	ARGS: 
		vel: float representing the current/initial velocity
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		atom_i: int representing the index of the current atom
		axis: int representing the axis to act on (0 - z) (1 - y) (2 - x)
		timestep: float
		mass: mass of the current atom
	RETURNS:
		float: projected acceleration on the atom

	This function calculates the expected acceleration of the atom over the
	given time. First, the displacement of the atom is calculated simply by
	using the input velocity, then the force is calculated using the 
	calc_single_force_RK4() method. The resulting force is then converted to 
	an acceleration in the appropriate units.
	''' 

	# Displacement of the current atom in specified coordinate axis
	disp = vel*timestep

	# Estimated geometry after timestep
	est_geometry = list()
	if(axis == 0):
		est_geometry = displace_geometry(geometry, disp,     0,    0, atom_i)
	elif(axis == 1):
		est_geometry = displace_geometry(geometry,     0, disp,    0, atom_i)
	elif(axis == 2):
		est_geometry = displace_geometry(geometry,     0,    0, disp, atom_i)
	else:
		sys.exit("\n\nAxis label not understood\n\n")

	force = calc_single_force_RK4(est_geometry, atom_i, axis)

	# F = ma   --->    a = F/m; also converts from m to angstrom
	return (force/mass)*10**10


def update_RK4(geometry, mass, atom_i, axis, velocities):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		atom_i: int representing the index of the current atom
		axis: int representing the axis to act on (0 - z) (1 - y) (2 - x)
		velocites: a list of floats representing the velocity of each atom 
		mass: mass of the current atom
	RETURNS:
		coord: coordinate of the atom in the current axis-direction (float)
		vel_f: final velocity of the atom in current direction (float)
	
	This method for RK4 will use the RK4 method to estimate the
		VELOCITY of the atom after time=dt. This essentially means that the 4
		terms calculated will be estimates of the ACCELERATION at various 
		points in the timestep. Then the velocity will be used to calculate
		the position by evolving the position as if the velocity was held 
		constant during this stretch. This is because the calculated velocity
		will be a weighted average for this timestep, already taking into 
		account its change over time. 
	''' 
	global dt
	# remember, geometry is a list of objects arranged as:
	# 							([atom] , ([z], [y], [x]))
	# so this is the current intial position of the atom in the appropriate axis
	coord = geometry[atom_i][1][axis]

	# velocities are stored in a 1D array and stored in order of z, y, x
	vel_i = velocities[(atom_i*3)+axis]

	v1 = vel_i
	v2 = vel_i + .5*calc_accel_RK4(v1, geometry, atom_i, axis, 0, mass)*dt
	v3 = vel_i + .5*calc_accel_RK4(v2, geometry, atom_i, axis, .5*dt, mass)*dt
	v4 = vel_i + .5*calc_accel_RK4(v3, geometry, atom_i, axis, .5*dt, mass)*dt

	# Estimated acceleration at the beginning of the interval
	a1 = calc_accel_RK4(v1, geometry, atom_i, axis, 0, mass)
	# First estimated acceleration at (dt/2)
	a2 = 2*calc_accel_RK4(v2, geometry, atom_i, axis, .5*dt, mass)
	# Second estimated acceleration at (dt/2)
	a3 = 2*calc_accel_RK4(v3, geometry, atom_i, axis, .5*dt, mass)
	# Estimate of acceleration at the end of the interval - dt
	a4 = calc_accel_RK4(v4, geometry, atom_i, axis, dt, mass)
	
	# Estimated velocity over interval by calculating the weighted average of
	# estimated accelerations using the formula tradtionally used in the 
	# Runge-Kutta 4th order method 
	vel_f = vel_i + dt*(a1+a2+a3+a4)/6

	coord += dt*vel_f

	return coord, vel_f


def runge_kutta_4(geometry, velocities, hamil=False):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		velocities: list of velocities of each atom in 3 axes (array of floats)
	RETURNS:
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		velocities: list of floats representing the updated velocities after
			the timestep

	This function uses the Runge-Kutta 4th order method to update the positions
	and velocities of each atom after a time of dt.
	''' 
	updated_locations = list()
	updated_velocities = list()
	for atom in range(N):
		# for each atom, find the updated coordinate and velocity in each axis
		mass = MASS_DICT[geometry[atom][0]]
		z_coord, z_vel = update_RK4(geometry, mass, atom, 0, velocities)
		y_coord, y_vel = update_RK4(geometry, mass, atom, 1, velocities)
		x_coord, x_vel = update_RK4(geometry, mass, atom, 2, velocities)
		# append the new geometry of the atom
		updated_locations.append((geometry[atom][0], 
			                      (z_coord, y_coord, x_coord)))
		updated_velocities.append(z_vel)
		updated_velocities.append(y_vel)
		updated_velocities.append(x_vel)
	if(hamil):
		ground_state_energy, hamiltonian = getGroundState(geometry)
		write_hamiltonians_to_file(hamiltonian)
	return updated_locations, updated_velocities


def euler_cromer(geometry, velocities, hamil=False):
	''' 
	ARGS: 
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		velocities: list of velocities of each atom in 3 axes (array of floats)
		hamil: boolean indicating whether to write the hamiltonian to disk
	RETURNS:
		geometry: list of tuples in format: 
			([atomic symbol], ([z], [y], [x]))
		velocities: list of floats representing the updated velocities after
			the timestep

	This function uses the Euler-Cromer method to update the positions
	and velocities of each atom after a time of dt.
	''' 
	# find the new forces given the current positions
	# average_forces should be a list of size N*3 (three forces per atom)
	average_forces = findForces_EC(geometry, hamil)

	# determine the new velocities by the forces acting on each atom
	velocities = update_velocities_EC(velocities, average_forces)

	# evolve the location of each atom to determine the geometry 
	geometry = update_positions_EC(geometry, velocities)

	return geometry, velocities


def evolve():
	''' 
	ARGS: 
		None
	RETURNS:
		None

	This function is the outermost loop of the simulation. Using the global
	variables INITIAL_GEOMETRY and INITIAL_VELOCITIES, set by the parse_input()
	function, this method loops over the number of iterations and writes the 
	new geometry to the file specified on input. Either the Euler-Cromer method
	or the Runge-Kutta 4th order method can be used in the inner loop and the 
	hamiltonian can be written to disk as well by using the modulo operator and 
	setting hamil=True for the desired method (RK4 implementation does not yet
	support this feature).
	''' 
	geometry = INITIAL_GEOMETRY
	velocities = INITIAL_VELOCITIES

	# write the current locations of each atom to the data file
	write_data(geometry)

	for x in range(0,ITERATIONS):
		# use the Euler-Cromer evolution method to approximate the next location
		# 		for each atom
		# if(x%10 == 0):
		# 	geometry, velocities = euler_cromer(geometry, velocities, 
		# 										hamil=True)
		# else:
		# 	geometry, velocities = euler_cromer(geometry, velocities)
		# if(x%10 == 0):
		# 	geometry, velocities = runge_kutta_4(geometry, velocities, hamil=True)
		# else:
		# 	geometry, velocities = runge_kutta_4(geometry, velocities)
		start = timeit.timeit()
		# geometry, velocities = runge_kutta_4(geometry, velocities)
		print(timeit.timeit()-start)
		# write the current locations of each atom to the data file
		write_data(geometry)


if __name__ == '__main__':
	parse_inputs(sys.argv)
	evolve()
