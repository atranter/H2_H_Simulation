from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

N_ATOMS = 0
N_ITERATIONS = 0

# Each object in the list is an atom. Each atom object has 3 lists -> 
	# z, y, x coords
ATOM_COORDS = list()


def get_data(filename):
	global N_ITERATIONS, N_ATOMS

	file = open(filename, "r")

	for line in file:

		N_ITERATIONS += 1

		coordinates = line.split(" ")
		N_ATOMS = len(coordinates)/3

		for i in range(N_ATOMS):
			ATOM_COORDS.append([list(), list(), list()])

		for i, coord in enumerate(coordinates):
			atom_i = i/3 # 3 coordinates per atom
			if(coord == ""):
				continue 
			if(i%3 == 0):
				ATOM_COORDS[atom_i][0].append(float(coord)) # -> atom's list
			elif(i%3 == 1):
				ATOM_COORDS[atom_i][1].append(float(coord))
			else:
				ATOM_COORDS[atom_i][2].append(float(coord))
	file.close()


# second option - move the point position at every frame
def update_point(n, x, y, z, point):
    point.set_data(np.array([x[n], y[n]]))
    point.set_3d_properties(z[n], 'z')
    return point


def plot():
	global x1, y1, z1, x2, y2, z2
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	for atom in ATOM_COORDS:
		atom[0] = np.asarray(atom[0])
		atom[1] = np.asarray(atom[1])
		atom[2] = np.asarray(atom[2])

	# create the first plot
	points = list()
	for atom in range(N_ATOMS):
		point, = ax.plot([ATOM_COORDS[atom][2][0]], [ATOM_COORDS[atom][1][0]], [ATOM_COORDS[atom][0][0]], 'o')
		points.append(point)

	ax.legend()
	ax.set_xlim([-2, 2])
	ax.set_ylim([-2, 2])
	ax.set_zlim([-2, 2])

	animations = list()
	for atom in range(N_ATOMS):
		animations.append(animation.FuncAnimation(fig, update_point, N_ITERATIONS, fargs=(ATOM_COORDS[atom][2], ATOM_COORDS[atom][1], ATOM_COORDS[atom][0], points[atom])))

	plt.show()


get_data("better.txt")
plot()
# show_data()
