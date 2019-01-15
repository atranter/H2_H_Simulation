from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation as ani
import sys

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


def update_lines(i, lines):
	global ATOM_COORDS, N_ATOMS
	j = 0
	if(i > 500):
		j = i-500
	for line, atom in zip(lines, range(N_ATOMS)):
		line.set_data([
					  ATOM_COORDS[atom][2][j:i],
					  ATOM_COORDS[atom][1][j:i]])
		line.set_3d_properties(ATOM_COORDS[atom][0][j:i])
	return lines


def plot(file):
	global N_ATOMS, N_ITERATIONS, ATOM_COORDS
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	lines = [ax.plot([ATOM_COORDS[atom][2][0]], 
					 [ATOM_COORDS[atom][1][0]], 
					 [ATOM_COORDS[atom][0][0]], '-', lw=1)[0] for atom in range(N_ATOMS)]

	ax.set_xlim([-1, 1])
	ax.set_ylim([-1, 1])
	ax.set_zlim([-1, 1])

	# Set up formatting for the movie files
	Writer = ani.writers['ffmpeg']
	writer = Writer(fps=100, metadata=dict(artist='Me'), bitrate=1800)

	line_animation = ani.FuncAnimation(
		fig, update_lines, N_ITERATIONS, fargs=([lines]), interval=1, 
		blit=True, repeat=True)

	file = str(file).replace('.txt', '.mp4')
	line_animation.save(file, writer=writer)
	plt.show()









if __name__ == '__main__':
	get_data(sys.argv[1])
	plot(sys.argv[1])