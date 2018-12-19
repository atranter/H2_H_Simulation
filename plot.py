from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

z1 = list()
y1 = list()
x1 = list()
z2 = list()
y2 = list()
x2 = list()
z3 = list()
y3 = list()
x3 = list()

def get_data(filename):
	file = open(filename, "r")
	for line in file:
		coordinates = line.split(" ")
		for i in range(0,9):
			if(i/3 == 0):
				if(i%3 == 0):
					z1.append(coordinates[i])
				elif(i%3 == 1):
					y1.append(coordinates[i])
				elif(i%3 == 2):
					x1.append(coordinates[i])
			elif(i/3 == 1):
				if(i%3 == 0):
					z2.append(coordinates[i])
				elif(i%3 == 1):
					y2.append(coordinates[i])
				elif(i%3 == 2):
					x2.append(coordinates[i])
			elif(i/3 == 2):
				if(i%3 == 0):
					z3.append(coordinates[i])
				elif(i%3 == 1):
					y3.append(coordinates[i])
				elif(i%3 == 2):
					x3.append(coordinates[i])

def show_data():
	for i in y2:
		print i

# second option - move the point position at every frame
def update_point(n, x, y, z, point):
    point.set_data(np.array([x[n], y[n]]))
    point.set_3d_properties(z[n], 'z')
    return point


def plot():
	global x1, y1, z1, x2, y2, z2, x3, y3, z3
	fig = plt.figure()
	ax = p3.Axes3D(fig)

	x1 = map(float,x1)
	y1 = map(float,y1)
	z1 = map(float,z1)

	x2 = map(float,x2)
	y2 = map(float,y2)
	z2 = map(float,z2)

	x3 = map(float,x3)
	y3 = map(float,y3)
	z3 = map(float,z3)

	x1 = np.asarray(x1)
	y1 = np.asarray(y1)
	z1 = np.asarray(z1)
	x2 = np.asarray(x2)
	y2 = np.asarray(y2)
	z2 = np.asarray(z2)
	x3 = np.asarray(x3)
	y3 = np.asarray(y3)
	z3 = np.asarray(z3)

	# create the first plot
	point1, = ax.plot([x1[0]], [y1[0]], [z1[0]], 'o')
	point2, = ax.plot([x2[0]], [y2[0]], [z2[0]], 'o')
	point3, = ax.plot([x3[0]], [y3[0]], [z3[0]], 'o')


	ax.legend()
	ax.set_xlim([-2, 2])
	ax.set_ylim([-2, 2])
	ax.set_zlim([-2, 2])

	ani1=animation.FuncAnimation(fig, update_point, 61, fargs=(x1, y1, z1, point1))
	ani2=animation.FuncAnimation(fig, update_point, 61, fargs=(x2, y2, z2, point2))
	ani3=animation.FuncAnimation(fig, update_point, 61, fargs=(x3, y3, z3, point3))

	plt.show()

get_data("better.txt")
plot()
# show_data()
