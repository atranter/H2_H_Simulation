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
		Variational Timestep
		Update to Z dimension for force finding'''

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


def write_hamiltonians_to_file(hamiltonian):
	filename = "hamil_" + DATAFILE
	file = open(filename, "a")
	file.write(str(hamiltonian))
	file.write("\n")
	file.write("\n")
	file.close()


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
    molecule.perform_transform("BK")
    # Convert from Hartree to eV
    # return 27.2114*molecule.molecule.hf_energy
    return 27.2114*molecule.ground_state_energy, molecule.qubit_hamiltonian_bk


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


def update_positions_EC(geometry, velocities):
	newGeometry = list()
	for i in range(N): # for each atom in our molecule, update the position
		z_coord = float(geometry[i][1][0]) + velocities[(i*3)]*dt
		y_coord = float(geometry[i][1][1]) + velocities[(i*3)+1]*dt
		x_coord = float(geometry[i][1][2]) + velocities[(i*3)+2]*dt

		newGeometry.append(tuple((str(geometry[i][0]), # atomic symbol
					             (z_coord, y_coord, x_coord)))) # position
	return newGeometry


def findForces_EC(geometry, hamil):
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
		energy, hamiltonian = getGroundState(g)
		energies.append(energy)
		if(hamil):
			write_hamiltonians_to_file(hamiltonian)

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


def update_velocities_EC(velocities, forces):
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


def calc_single_force_RK4(geometry, atom_i, axis):
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


def calc_velocity_RK4(timestep, current_vel, force, mass):
	return current_vel + (force/mass)*timestep*10**10


def update_RK4(geometry, mass, atom_i, axis, velocities):
	global dt
	coord = geometry[atom_i][1][axis]
	current_vel = velocities[(atom_i*3)+axis]
	force = calc_single_force_RK4(geometry, atom_i, axis)

	# k1v will always be the current_vel since there will be no more accel
	k1v = calc_velocity_RK4(0,          current_vel,              force, mass)
	# k2v = current_vel + (force/mass)*(dt/2)*10**10
	k2v = calc_velocity_RK4(0 + (dt/2), current_vel + (k1v*dt/2), force, mass)
	k3v = calc_velocity_RK4(0 + (dt/2), current_vel + (k2v*dt/2), force, mass)
	k4v = calc_velocity_RK4(0 + dt,     current_vel + (k3v*dt),   force, mass)

	velocity = current_vel + (k1v + 2*k2v + 2*k3v + k4v)/6

	# k1 = current_vel
	# k2 = current_vel + k1v*dt/2
	# k3 = current_vel + k2v*dt/2
	# k4 = current_vel + k3v*dt

	# coord += (k1 + 2*k2 + 2*k3 + k4)*(dt/6)

	coord += dt*velocity

	# NOTE, MUST STORE CURRENT VELOCITIES
	# k1 = dt*velocity(0, coord, force)
	# k2 = dt*velocity((0.5*dt), coord+(0.5*k1), force, mass, current_vel)
	# k3 = dt*velocity((0.5*dt), coord+(0.5*k2), force, mass, current_vel)
	# k4 = dt*velocity(dt,       coord+k3,       force, mass, current_vel)

	# coord += (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)

	return coord, velocity


def calc_accel_RK4(vel, geometry, atom_i, axis, timestep, mass):
	''' In order to estimate the acceleration at the given time, we need to:
			1. Have the velocity, geometry, atom_i, axis, timestep, mass
			2. Using the given timestep and intial velocity, estimate the 
				position of the target atom after timestep and update the geo
			3. Using this new geometry, calculate the average force on the atom
				in the given coordinate axis
			4. Return force/mass = acceleration '''
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



def new_update_RK4(geometry, mass, atom_i, axis, velocities):
	''' This method for RK4 will basically use the RK4 method to estimate the
		VELOCITY of the atom after time=dt. This essentially means that the 4
		terms calculated will be estimates of the ACCELERATION at various 
		points in the timestep. Then the velocity will be used to calculate
		the position by evolving the position as if the velocity was held 
		constant during this stretch. This is because the calculated velocity
		will be a weighted average for this timestep, already taking into 
		account its change over time. '''
	global dt
	# remember, geometry is a list of objects arranged as:
	# 							([atom] , ([z], [y], [x]))
	coord = geometry[atom_i][1][axis]
	# velocities are stored in a 1D array and stored in order of z, y, x
	vel_i = velocities[(atom_i*3)+axis]

	v1 = vel_i
	v2 = vel_i + .5*calc_accel_RK4(v1, geometry, atom_i, axis, 0, mass)*dt
	v3 = vel_i + .5*calc_accel_RK4(v2, geometry, atom_i, axis, .5*dt, mass)*dt
	v4 = vel_i + .5*calc_accel_RK4(v3, geometry, atom_i, axis, .5*dt, mass)*dt

	a1 = calc_accel_RK4(v1, geometry, atom_i, axis, 0, mass)
	a2 = 2*calc_accel_RK4(v2, geometry, atom_i, axis, .5*dt, mass)
	a3 = 2*calc_accel_RK4(v3, geometry, atom_i, axis, .5*dt, mass)
	a4 = calc_accel_RK4(v4, geometry, atom_i, axis, dt, mass)

	vel_f = vel_i + dt*(a1+a2+a3+a4)/6

	# # Estimated acceleration at the beginning of the interval
	# a1 = calc_accel_RK4(vel_i,           geometry, atom_i, axis,      0, mass)
	# # First estimated acceleration at (dt/2)
	# a2 = calc_accel_RK4(vel_i+(a1*dt/2), geometry, atom_i, axis, (dt/2), mass)
	# # Second estimated acceleration at (dt/2)
	# a3 = calc_accel_RK4(vel_i+(a2*dt/2), geometry, atom_i, axis, (dt/2), mass)
	# # Estimate of acceleration at the end of the interval - dt
	# a4 = calc_accel_RK4(vel_i+(a3*dt),   geometry, atom_i, axis,     dt, mass)

	# # Estimated velocity over interval by calculating the weighted average of
	# # estimated accelerations using the formula tradtionally used in the 
	# # Runge-Kutta 4th order method 
	# vel_f = dt*(a1+2*a2+2*a3+a4)/6

	print vel_f

	coord += dt*vel_f

	return coord, vel_f


def runge_kutta_4(geometry, velocities):
	updated_locations = list()
	updated_velocities = list()
	for atom in range(N):
		mass = MASS_DICT[geometry[atom][0]]
		# Make this multi-threaded
		z_coord, z_vel = new_update_RK4(geometry, mass, atom, 0, velocities)
		y_coord, y_vel = new_update_RK4(geometry, mass, atom, 1, velocities)
		x_coord, x_vel = new_update_RK4(geometry, mass, atom, 2, velocities)
		updated_locations.append((geometry[atom][0], 
			                      (z_coord, y_coord, x_coord)))
		updated_velocities.append(z_vel)
		updated_velocities.append(y_vel)
		updated_velocities.append(x_vel)
	geometry = updated_locations
	velocities = updated_velocities
	return geometry, velocities


def euler_cromer(geometry, velocities, hamil=False):
	# find the new forces given the current positions
	# average_forces should be a list of size N*3 (three forces per atom)
	average_forces = findForces_EC(geometry, hamil)

	# determine the new velocities by the forces acting on each atom
	velocities = update_velocities_EC(velocities, average_forces)

	# evolve the location of each atom to determine the geometry 
	geometry = update_positions_EC(geometry, velocities)

	return geometry, velocities


def evolve():
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


		# trying to make Runge-Kutta 4th Order method work
		geometry, velocities = runge_kutta_4(geometry, velocities)

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
	evolve()

if __name__ == '__main__':
	main(sys.argv)
