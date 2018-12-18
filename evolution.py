from OpenFermionWrapper import OpenFermionWrapper
import sys

INITIAL_VELOCITY_X = 0 #angstroms/s
INITIAL_VELOCITY_Y = 0 #angstroms/s
dt = 1000000000*10**-15 #s from fs
dx = 0.001 #angstrom
dy = 0.001 #angstrom
M = 1000000*938.2720813 #eV/c^2

INITIAL_GEOMETRY = [("H", (0,0,0)), ("H", (0,0,0.7414)), ("H", (0,0,2))]
BASIS = "sto-3g"
MULTIPLICITY = 2
CHARGE = 0

MIN_FORCE = 1000000000

positions = list()

def getGroundState(geometry):
	molecule = OpenFermionWrapper()
	molecule.load_molecule(geometry=geometry,         basis=BASIS, 
						   multiplicity=MULTIPLICITY, charge=CHARGE,
						   forceCalculation=True)
	molecule.set_ground_state_energy()
	# Convert from Hartree to eV
	return 27.2114*molecule.ground_state_energy

def newGeometry(geometry, dx, dy, i):
	newGeometry = list()
	if(i != 0):
		newGeometry.append(tuple((str(geometry[0][0]), 
					(float(geometry[0][1][0]), float(geometry[0][1][1]),
				     float(geometry[0][1][2])))))
	else:
		newGeometry.append(tuple((str(geometry[0][0]),
					(float(geometry[0][1][0]), float(geometry[0][1][1]+dy),
					 float(geometry[0][1][2]+dx)))))
	if(i != 1):
		newGeometry.append(tuple((str(geometry[1][0]), 
					(float(geometry[1][1][0]), float(geometry[1][1][1]),
					 float(geometry[1][1][2])))))
	else:
		newGeometry.append(tuple((str(geometry[0][0]),
					(float(geometry[1][1][0]), float(geometry[1][1][1]+dy),
					 float(geometry[1][1][2]+dx)))))
	if(i != 2):
		newGeometry.append(tuple((str(geometry[2][0]), 
					(float(geometry[2][1][0]), float(geometry[2][1][1]),
					 float(geometry[2][1][2])))))
	else:
		newGeometry.append(tuple((str(geometry[2][0]), 
					(float(geometry[2][1][0]), float(geometry[2][1][1]+dy),
					 float(geometry[2][1][2]+dx)))))
	return newGeometry

def findForces(geometry):
	geometries = list()
	energies = list()
	forces = list()
	average_forces = list()

	geometries.append(geometry)

	for i in range(0, 3):
		_r_geometry = newGeometry(geometry, dx, 0, i)
		geometries.append(_r_geometry)
		_l_geometry = newGeometry(geometry, -dx, 0, i)
		geometries.append(_l_geometry)
		_u_geometry = newGeometry(geometry, 0, dy, i)
		geometries.append(_u_geometry)
		_d_geometry = newGeometry(geometry, 0, -dy, i)
		geometries.append(_d_geometry)

	for g in geometries:
		# Energy given in eV
		energies.append(getGroundState(g))

	for e in energies:
		# Store as eV/m
		forces.append(-1*(energies[0]-e)/(dx*10**-10))

	for i in range(0, 3):
		# index is i*3+1?
		average_x = (forces[(4*i)+1]-forces[(4*i)+2])/2
		average_y = (forces[(4*i)+3]-forces[(4*i)+4])/2
		average_forces.append(average_x)
		average_forces.append(average_y)

	return average_forces

def evolveAll():
	geometry = INITIAL_GEOMETRY
	velocity_x = INITIAL_VELOCITY_X
	velocity_y = INITIAL_VELOCITY_Y
	velocities = list()
	for i in range(0, 4):
		velocities.append(0)
	velocities.append(velocity_x)
	velocities.append(velocity_y)

	file = open("evolutionall.txt", "w")

	average_forces = findForces(geometry)
	file.write("{} {} {} {} {} {} {} {} {}".format(average_forces[0], average_forces[1],
									   average_forces[2], average_forces[3],
									   average_forces[4], average_forces[5],
									   velocity_x, velocity_y, geometry)) 

	for x in range(0,15):
		for i in range(0, 3):
			velocity_x = velocities[(i*2)]
			velocity_y = velocities[(i*2)+1]
			fx = average_forces[(i*2)]
			fy = average_forces[(i*2)+1]
			velocity_x += (fx/M)*dt*10**10
			velocity_y += (fy/M)*dt*10**10
			velocities[(i*2)] = velocity_x
			velocities[(i*2)+1] = velocity_y

			geometry = newGeometry(geometry, velocity_x*dt, velocity_y*dt, i)

		average_forces = findForces(geometry)
		file.write("{} {} {} {} {} {} {} {} {} {}".format(average_forces[0], average_forces[1],
									   average_forces[2], average_forces[3],
									   average_forces[4], average_forces[5],
									   velocity_x, velocity_y, geometry)) 
	file.close()


def evolve():
	geometry = INITIAL_GEOMETRY
	positions.append(geometry)

	file = open("evolution.txt", "w")

	# forces are returned in eV/m
	average_forces = findForces(geometry)
	fx = average_forces[4]
	fy = average_forces[5]

	velocity_x, velocity_y = INITIAL_VELOCITY_X, INITIAL_VELOCITY_Y
	print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n ----- STATS INCOMING ---- \n")
	file.write("{} {} {} {} {}".format(fx, fy, velocity_x, velocity_y, geometry)) 
	print("\n\n ----- ENDING STATS ----- \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
	while(abs(fx) > MIN_FORCE or abs(fy) > MIN_FORCE):
		velocity_x += (fx/M)*dt*10**10 # converting to angstroms/s
		velocity_y += (fy/M)*dt*10**10

		geometry = newGeometry(geometry, velocity_x*dt, velocity_y*dt, 2)
		positions.append(geometry)

		# forces are returned in eV/m
		average_forces = findForces(geometry)
		fx = average_forces[4]
		fy = average_forces[5]
		print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n ----- STATS INCOMING ---- \n")
		file.write("{} {} {} {} {}".format(fx, fy, velocity_x, velocity_y, geometry)) 
		print("\n\n ----- ENDING STATS ----- \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")

	file.close()

def main(input):
	# evolve()
	evolveAll()
	# for position in positions:
	# 	print position
	# e = getGroundState([("H", (0,0,0))])
	# fx, fy = findForces(INITIAL_GEOMETRY)
	# print fx, fy
	# print 27.2114*e

if __name__ == '__main__':
	main(sys.argv)