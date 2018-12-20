import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import math

def generate_atoms():
	r = 3.822  #r_min for argon
	x_position = []
	y_position = []
	atom_position = []
	coor = 0.0
	for i in range(20):
		for j in range(20):
			x_position.append(coor)
		y_position.append(coor)
		coor += r
	y_position *= 20
	atom_position.append(x_position)
	atom_position.append(y_position)

	return atom_position
	#print(x_position)
    #r = sqrt((x(i,1)-x(j,1))**2 + (x(i,2)-x(j,2))**2)

def plot_atoms(atom_position):
	x = atom_position[0]
	y = atom_position[1]
	plt.plot(x,y,'ro')
	plt.xlabel("x_position (Å)")
	plt.ylabel("y_position (Å)")
	plt.show()

def neighbor_list(atom_position):
	'''atom_position is the 
	(2,400) array that contains the positions of
	the 400 atoms block
	This function returns the list of neighbors map
	of each atom in the block
    And a numNeigh list that contains the number 
    of neighbor of each atom
	'''

	cut_off = 7.5
	neighbors = []
	numNeigh = []
	for i in range(400):
		lst = []
		for j in range(400):
			if j != i:
				r = math.sqrt((atom_position[0][i] - atom_position[0][j])**2 + (atom_position[1][i]-atom_position[1][j])**2)
				if r < cut_off:
					lst.append(j)
		neighbors.append(lst)
	for i in range(len(neighbors)):
		numNeigh.append(len(neighbors[i]))
	
	#print(atom_position[0][21])
	#print(atom_position[1][21])
	#print(atomic_position)
	#print(numNeigh)
	return neighbors, numNeigh
	

def plot_neighbors(atom_position):
	x = atom_position[0]
	y = atom_position[1]
	fig, ax = plt.subplots()
	atoms = plt.plot(x,y,'ro')
	plt.xlabel("x_position (Å)")
	plt.ylabel("y_position (Å)")
	circle = plt.Circle((atom_position[0][300],atom_position[1][30]),7.5,color = 'blue')
	ax.add_artist(circle)
	plt.show()

def calc_atomic_level_stresses(atom_position,neighbors):
	sigma = 3.405
	epsilon = 0.010323
	x = atom_position[0]
	y = atom_position[1]
	F_X = []
	F_Y = []
	energy = []
	e = 0
	Force = []
	for i in range(400):
		ex = 0
		ey = 0
		for j in range(len(neighbors[i])):
			r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
			if r == 0:
				e = e
			else:
				e += 4*epsilon*((sigma/r)**12-(sigma/r)**6)
			 
		energy.append(e)
	r = np.arange(0, 72.618,72.618/len(energy))
	plt.plot(r,energy)
	plt.title("Energy vs Radius")
	plt.xlabel("Radius (Å)")
	plt.ylabel("Energy (eV)")
	plt.show()

	for i in range(400):
		fx = 0
		fy = 0
		for j in range(len(neighbors[i])):
			r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
			if x[neighbors[i][j]]-x[i] == 0:
				fx +=0
			else:
				fx += -24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
			if y[neighbors[i][j]]-y[i] == 0:
				fy +=0
			else:
				fy += -24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
		F_X.append(fx)
		F_Y.append(fy)
	Force.append(F_X)
	Force.append(F_Y)
	#print(Force)
	

	#Voronoi polyhedra
	volume = []
	for i in range(400):
		if x[i] == min(x) or y[i] == min(x) \
		or x[i] == max(x) or y[i] == max(y):
			v = 1.911*3.822
		if x[i] == min(x) and y[i] == min(y):
			v = 1.911**2
		if x[i] == min(x) and y[i] == min(y):
			v = 1.911**2
		if x[i] == max(x) and y[i] == min(y):
			v = 1.911**2
		if x[i] == max(x) and y[i] == max(y):
			v = 1.911**2
		elif x[i] != min(x) and x[i] != max(x) \
		and y[i] != max(y) and y[i] != min(y):
			v = 3.822**2
		volume.append(v)

	sigma11 = [] #xx
	sigma12 = [] #xy
	sigma22 = [] #yy
	for i in range(400):
		sum11 = 0
		sum12 = 0
		sum22 = 0
		for j in range(len(neighbors[i])):
			rx = atom_position[0][neighbors[i][j]] - atom_position[0][i]
			ry = atom_position[1][neighbors[i][j]] - atom_position[1][i]
			fx = Force[0][i]
			fy = Force[1][i]
			sum11 += fx*rx
			sum12 += fx*ry
			sum22 += fy*ry

		sig11 = -1/volume[i]*sum11
		sig12 = -1/volume[i]*sum12
		sig22 = -1/volume[i]*sum22
		sigma11.append(sig11)
		sigma12.append(sig12)
		sigma22.append(sig22)


	fig = plt.scatter(x,y,c=sigma11) # plot with coloring by the first column
	a = plt.colorbar() # create colorbar
	a.set_label("${\sigma_{11}}$ (eV/${\AA^2}$)") # label colorbar
	#plt.savefig("sigma11.pdf") # write to file
	plt.title('${\sigma_{11}}$')
	plt.xlabel('x (Å)')
	plt.ylabel('y (Å)')
	plt.show() # display the plot

	fig = plt.scatter(x,y,c=sigma12) # plot with coloring by the first column
	a = plt.colorbar() # create colorbar
	a.set_label("${\sigma_{12}}$ (eV/${\AA^2}$)") # label colorbar
	plt.title('${\sigma_{12}}$')
	plt.xlabel('x (Å)')
	plt.ylabel('y (Å)')
	plt.show() # display the plot

	fig = plt.scatter(x,y,c=sigma22) # plot with coloring by the first column
	a = plt.colorbar() # create colorbar
	a.set_label("${\sigma_{22}}$ (eV/${\AA^2}$)") # label colorbar
	plt.title('${\sigma_{22}}$')
	plt.xlabel('x (Å)')
	plt.ylabel('y (Å)')
	plt.show() # display the plot




def calc_RDF(delta_r,atom_position,neighbors):
	V = 72.618**2
	N = 400
	pho = N/V
	pi = 3.1415926
	r = 0.5
	gr = []
	lst = []
	while r < 15:
		count = 0
		nr = 0
		for i in range(400):
			for j in range(400):
				radius = math.sqrt((atom_position[0][i] - atom_position[0][j])**2 + (atom_position[1][i]-atom_position[1][j])**2)
				if r < radius < r + delta_r:
					count += 1
					nr += 1/N*count
		gr.append(1/pho*nr/(2*pi*r*delta_r))
		lst.append(r)
		r += 0.01

	plt.plot(lst,gr)
	plt.title("RDF vs Radius for delta_r = %s" %(delta_r))
	plt.xlabel("Radius (Å)")
	plt.ylabel("RDF g(r)")
	plt.show()


	


def main():
	atom_position = generate_atoms()
	plot_atoms(atom_position)
	neighbors, numNeigh = neighbor_list(atom_position)
	plot_neighbors(atom_position)
	calc_atomic_level_stresses(atom_position,neighbors)
	calc_RDF(0.01,atom_position,neighbors)

if __name__ == '__main__':
	main()