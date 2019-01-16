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


def plot_lines(file, save=False):
	global N_ATOMS, N_ITERATIONS, ATOM_COORDS
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	lines = [ax.plot([ATOM_COORDS[atom][2][0]], 
					 [ATOM_COORDS[atom][1][0]], 
					 [ATOM_COORDS[atom][0][0]], '-', lw=1)[0] for atom in range(N_ATOMS)]

	ax.set_xlim([-1, 1])
	ax.set_ylim([-1, 1])
	ax.set_zlim([-1, 1])

	line_animation = ani.FuncAnimation(
		fig, update_lines, N_ITERATIONS, fargs=([lines]), interval=1, 
		blit=True, repeat=True)

	if(save):
		# Set up formatting for the movie files
		Writer = ani.writers['ffmpeg']
		writer = Writer(fps=100, metadata=dict(artist='Me'), bitrate=1800)
		file = str(file).replace('.txt', '_lines.mp4')
		line_animation.save(file, writer=writer)

	plt.show()


def update_points(i, points):
	global ATOM_COORDS, N_ATOMS
	for point, atom in zip(points, range(N_ATOMS)):
		point.set_data([
					  ATOM_COORDS[atom][2][i],
					  ATOM_COORDS[atom][1][i]])
		point.set_3d_properties(ATOM_COORDS[atom][0][i])
	return points


def plot_points(file, save=False):
	global N_ATOMS, N_ITERATIONS, ATOM_COORDS
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	points = [ax.plot([], [], [], 'o')[0] for atom in range(N_ATOMS)]

	ax.set_xlim([-1, 1])
	ax.set_ylim([-1, 1])
	ax.set_zlim([-1, 1])

	point_animation = ani.FuncAnimation(
		fig, update_points, N_ITERATIONS, fargs=([points]), interval=1, 
		blit=True, repeat=True)

	if(save):
		# Set up formatting for the movie files
		Writer = ani.writers['ffmpeg']
		writer = Writer(fps=100, metadata=dict(artist='Me'), bitrate=1800)
		file = str(file).replace('.txt', '_points.mp4')
		line_animation.save(file, writer=writer)
	
	plt.show()


if __name__ == '__main__':
	save = False
	if(len(sys.argv) == 3):
		if(sys.argv[2] == "y"):
			save = True
		elif(sys.argv[2] == "n"):
			save = False
		else:
			sys.exit('''\n\n----Error: Could not understand save command: {} ----\n\n'''
			   .format(sys.argv[2]))
	elif(len(sys.argv) != 2):
		sys.exit("\n\n----Error: Incorrect number of input parameters----\n\n")

	get_data(sys.argv[1])
	plot_lines(sys.argv[1], save=save)
	# plot_points(sys.argv[1], save=save)