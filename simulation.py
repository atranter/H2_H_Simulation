''' INPUT FORMAT (all on one line):
			   python simulation.py [datafilename] [multiplicity] [charge]/
for each atom: [atomic symbol] [z coord] [y coord] [x coord]/
			   [velocity z] [velocity y] [velocity x]

	EXAMPLE FORMAT:
	python simulation.py example.txt 1 0 H 0 0 0 0 0 0 H 0 0 0.7414 0 0 0



	TODO:
		Runge-Kutta 4
			Q: do we set dxdt aka velocity to be directly determined from force?
			or do we recalculate force everytime?
		Variational Timestep'''

from OpenFermionWrapper import OpenFermionWrapper
import sys
import math

''' --- Set Constants ---'''
dt = 10000000*10**-15 #s from fs
dx = 0.001 #angstrom
dy = 0.001 #angstrom -> how much a molecule is displaced to find nearby energies
dz = 0.001 #angstrom
M = 1000000*938.2720813 #eV/c^2 -> this needs to become a dictionary
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
DATAFILE = "uhoh"
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

	if(((len(input)-4)%7 != 0) or (len(input) < 11)):
		sys.exit("Incorrect number of input parameters. See documention.\n\n")

	DATAFILE = input[1]
	MULTIPLICITY = int(input[2])
	CHARGE = int(input[3])
	N = ((len(input)-4)/7)

	for i in range(N):
		atom_i = (i*7)+4
		# print MASS_DICT[str(input[atom_i])]
		MASSES.append(MASS_DICT[str(input[atom_i])])
		atom = tuple((str(input[atom_i]), (float(input[atom_i+1]), 
				   float(input[atom_i+2]), float(input[atom_i+3]))))
		INITIAL_GEOMETRY.append(atom)
		INITIAL_VELOCITIES.append(float(input[atom_i+4]))
		INITIAL_VELOCITIES.append(float(input[atom_i+5]))
		INITIAL_VELOCITIES.append(float(input[atom_i+6]))


def getGroundState(geometry):
    ''' 
    ARGS: 
    
    This function obtains the ground state energy of the given molecule
    	using the OpenFermionWrapper class, the OpenFermion program, and
    	the Psi4 program. It then converts from units of Hartree and returns
    	the energy of the molecule in eV.'''
    molecule = OpenFermionWrapper()
    molecule.load_molecule(geometry=geometry,         basis=BASIS, 
    					   multiplicity=MULTIPLICITY, charge=CHARGE,
    					   forceCalculation=True)
    molecule.set_ground_state_energy()
    # Convert from Hartree to eV
    return 27.2114*molecule.molecule.hf_energy
    # 	return 27.2114*molecule.ground_state_energy


def update_positions_EC(geometry, velocities):
	newGeometry = list()
	for i in range(N): # for each atom in our molecule, update the position
		z_coord = float(geometry[i][1][0]) + velocities[(i*3)]*dt
		y_coord = float(geometry[i][1][1]) + velocities[(i*3)+1]*dt
		x_coord = float(geometry[i][1][2]) + velocities[(i*3)+2]*dt

		newGeometry.append(tuple((str(geometry[i][0]), # atomic symbol
					             (z_coord, y_coord, x_coord)))) # position
	return newGeometry


def displace_geometry(geometry, dz, dy, dx, a):
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


def findForces(geometry):
	geometries = list()
	energies = list()
	forces = list()
	average_forces = list()

	geometries.append(geometry)

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

	for g in geometries:
		# Energy given in eV
		energies.append(getGroundState(g))

	for e in energies:
		# Store as eV/m
		forces.append(-1*(e-energies[0])/(dx*10**-10))

	for i in range(N):
		# index is i*3+1?
		average_x = (forces[(4*i)+1]-forces[(4*i)+2])/2
		average_y = (forces[(4*i)+3]-forces[(4*i)+4])/2
		average_forces.append(0) # not worrying about forces in z yet
		average_forces.append(average_y)
		average_forces.append(average_x)

	return average_forces


def write_data(geometry):
	data = open(DATAFILE, "a")
	for i in range(N):
		if(i == N-1):
			data.write("{} {} {}".format(geometry[i][1][0], geometry[i][1][1],
									     geometry[i][1][2]))
		else:
			data.write("{} {} {} ".format(geometry[i][1][0], geometry[i][1][1],
									      geometry[i][1][2]))
	data.write("\n")
	data.close()
	
	
def get_time_step(vi, d, a):
    vf = math.sqrt((vi**2) + 2*a*d)
    time = (2*d)/(vi+vf)
    return time


def set_time_step(velocities, forces):
    global dt
    index = forces.index(max(forces, key=abs))
    mass = MASSES[index/3]
    a = (forces[index]/mass)*10**10
    dt = get_time_step(velocities[index], 0.01, a)


def update_velocities(velocities, forces):
    # set_time_step(velocities, forces)
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


# NEEDS: timestep, geometry, coordinate index, atom index
def velocity(timestep, geometry, axis, atom):
	# FUCK... have to find the force from here
	# STEP 1: using geometry, find force
	# STEP 2: find current velocity from velocities array
	# STEP 3: set current velocity in velocities array to new velocity
	#         calculated after time=dt
	# STEP 4: return the estimated velocity after time=timestep
	return current_vel + (force/mass)*timestep**10**10


def runge_kutta_4(coord, force, mass):
	global dt
	# NOTE, MUST STORE CURRENT VELOCITIES
	k1 = dt*velocity(0, coord, force)
	k2 = dt*velocity((0.5*dt), coord+(0.5*k1), force, mass, current_vel)
	k3 = dt*velocity((0.5*dt), coord+(0.5*k2), force, mass, current_vel)
	k4 = dt*velocity(dt,       coord+k3,       force, mass, current_vel)

	coord += (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

	return coord


def update_positions_RK(geometry, time, forces):
	# time is given as the xth iteration, so it can then be calculated using dt
	global dt
	time *= dt
	updated_locations = list()
	for atom in range(N):
		mass = MASS_DICT[geometry[atom][0]]
		z_coord = runge_kutta_4(time, geometry[atom][1][0], h, mass)
		y_coord = runge_kutta_4(time, geometry[atom][1][1], h, mass)
		x_coord = runge_kutta_4(time, geometry[atom][1][2], h, mass)
		updated_locations.append(geometry[atom][0], (z_coord, y_coord, x_coord))
	return updated_locations


def evolveAll():
	geometry = INITIAL_GEOMETRY
	velocities = INITIAL_VELOCITIES

	# write the current locations of each atom to the data file
	write_data(geometry)

	for x in range(0,ITERATIONS):
		# find the new forces given the current positions
		# average_forces should be a list of size N*3 (three forces per atom)
		average_forces = findForces(geometry)
	
		# determine the new velocities by the forces acting on each atom
		velocities = update_velocities(velocities, average_forces)

		# evolve the location of each atom to determine the geometry 
		geometry = update_positions_EC(geometry, velocities)

		# write the current locations of each atom to the data file
		write_data(geometry)


def check_parse():
	print DATAFILE
	print MULTIPLICITY
	print CHARGE
	print N
	for atom in INITIAL_GEOMETRY:
		print atom
	for vel in INITIAL_VELOCITIES:
		print vel


def main(input):
	parse_inputs(input)
	# check_parse()
	evolveAll()

if __name__ == '__main__':
	main(sys.argv)