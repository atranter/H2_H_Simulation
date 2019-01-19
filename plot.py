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
	''' 
	ARGS:
		filename - a string with the local location of the file with the stored
				   data
	RETURNS:
		None (sets the global variables N_ATOMS, N_ITERATIONS, ATOM_COORDS)

	This function parses the input file into the data structures used in this
	script. The datafile must contain the positions of atoms in the format
		[atom1_z] [atom1_y] [atom1_x] ... for every atom in the molecule
	Each new line is the updated coordinate positions after the next timestep.

	'''
	global N_ITERATIONS, N_ATOMS

	file = open(filename, "r")

	for line in file:
		# each line is the geometry at a new point in time 
		N_ITERATIONS += 1

		coordinates = line.split(" ")

		if(N_ATOMS == 0):
			# each atom has z, y, x coordinates, so the number of atoms follows from
			# the number of coordinates on the line
			N_ATOMS = len(coordinates)/3

		if(ATOM_COORDS == []):
			# if the ATOM_COORDS list is empty, add an object with 3 empty lists
			# for each atom to store the z, y, x coordinates at each iteration
			for i in range(N_ATOMS):
				ATOM_COORDS.append([list(), list(), list()])

		for i, coord in enumerate(coordinates):
			atom_i = i/3 # 3 coordinates per atom
			if(coord == ""):
				continue 
			# append the coordinate to the appropriate atom and axis list after
			# casting it as a float
			if(i%3 == 0):
				ATOM_COORDS[atom_i][0].append(float(coord))
			elif(i%3 == 1):
				ATOM_COORDS[atom_i][1].append(float(coord))
			else:
				ATOM_COORDS[atom_i][2].append(float(coord))
	file.close()


def update_lines(i, lines):
	''' 
	ARGS:
		i - integer index of the current iteration
		lines - an array of ax.plot() objects that represent lines
	RETURNS:
		lines - updated array of ax.plot() objects
	
	This function updates the lines in the animation to the next points. After
	500 iterations, the end of the line will begin to disappear as the front
	continues to move (only a maximum of 500 points are included in the line). 
	'''
	global ATOM_COORDS, N_ATOMS

	# j variable keeps track of the beginning of the line
	j = 0
	if(i > 500):
		# If you want to make the line shorter or larger, just decrease or
		# increase from 500 respectively
		j = i-500

	for line, atom in zip(lines, range(N_ATOMS)):
		# update the data that goes into each line
		line.set_data([
					  ATOM_COORDS[atom][2][j:i],
					  ATOM_COORDS[atom][1][j:i]])
		line.set_3d_properties(ATOM_COORDS[atom][0][j:i])
	return lines


def plot_lines(file, save=False):
	'''
	ARGS: 
		file - string with the local location of the data file. Will be used 
			   to save the movie in mp4 format
		save - a boolean (preset to False) that determines whether or not to
			   save the movie
	RETURNS:
		none

	This function plots the evolution of the molecule over time using lines. 
	A line for every atom is used to more clearly show the path of the atoms
	over time. A plot of this animation is shown through the display (client
	must make sure to set the display variables correctly) and is saved to mp4
	format depending on the value of the save variable. This function primarily
	depends on the matplotlib plots and animations. 
	'''
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
	''' 
	ARGS:
		i - integer index of the current iteration
		points - an array of ax.plot() objects that represent points
	RETURNS:
		points - updated array of ax.plot() objects
	
	This function updates the points in the animation to the next locations.
	'''
	global ATOM_COORDS, N_ATOMS
	for point, atom in zip(points, range(N_ATOMS)):
		point.set_data([
					  ATOM_COORDS[atom][2][i],
					  ATOM_COORDS[atom][1][i]])
		point.set_3d_properties(ATOM_COORDS[atom][0][i])
	return points


def plot_points(file, save=False):
	'''
	ARGS: 
		file - string with the local location of the data file. Will be used 
			   to save the movie in mp4 format
		save - a boolean (preset to False) that determines whether or not to
			   save the movie
	RETURNS:
		none

	This function plots the evolution of the molecule over time using points. 
	A point of uniform size is used for each atom. A plot of this animation is 
	shown through the display (client must make sure to set the display 
	variables correctly) and is saved to mp4 format depending on the value of 
	the save variable. This function primarily depends on the matplotlib plots 
	and animations. 
	'''
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
	# plot_lines(sys.argv[1], save=save)
	plot_points(sys.argv[1], save=save)