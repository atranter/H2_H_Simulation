from OpenFermionWrapper import OpenFermionWrapper
import sys

INITIAL_VELOCITY_Y = 0 #angstroms/s
dt = 20000000*10**-15 #s from fs
INITIAL_VELOCITY_X = 0/dt #angstroms/s
dx = 0.001 #angstrom
dy = 0.001 #angstrom
M = 1000000*938.2720813 #eV/c^2

INITIAL_GEOMETRY = [("H", (0,0,0)), ("H", (0,0,0.7414))]
N = len(INITIAL_GEOMETRY)
BASIS = "sto-3g"
MULTIPLICITY = 1
CHARGE = 0

ITERATIONS = 50

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
	return newGeometry

def findForces(geometry):
	geometries = list()
	energies = list()
	forces = list()
	average_forces = list()

	geometries.append(geometry)

	for i in range(0, 2):
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

	for i in range(0, 2):
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
	for i in range(0, 2):
		velocities.append(0)
	velocities.append(velocity_x)
	velocities.append(velocity_y)

	average_forces = findForces(geometry)

	data = open("evolution_H2Sim_data.txt", "w")
	data.write("{} {} {} {} {} {}\n".format(geometry[0][1][0], 
											geometry[0][1][1],
											geometry[0][1][2], 
											geometry[1][1][0],
											geometry[1][1][1], 
											geometry[1][1][2]))
	data.close()

	for x in range(0,ITERATIONS):
		for i in range(0, 2):
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
		data = open("evolution_H2Sim_data.txt", "a")
		data.write("{} {} {} {} {} {}\n".format(geometry[0][1][0], 
												geometry[0][1][1],
												geometry[0][1][2], 
												geometry[1][1][0],
												geometry[1][1][1], 
												geometry[1][1][2]))
		data.close()

def main(input):
	evolveAll()

if __name__ == '__main__':
	main(sys.argv)