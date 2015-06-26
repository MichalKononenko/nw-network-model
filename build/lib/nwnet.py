#! /usr/bin/env python
"""
nwnet.py
Random resitor network module

Created by Jeremy Smith on 2015-06-07
University of California, Berkeley
j-smith@ecs.berkeley.edu
"""

import time
import numpy as np
from scipy import linalg
import turtle as tu

__author__ = "Jeremy Smith"
__version__ = "2.2"


# Useful function definitions
def intersectCheck(start1, end1, start2, end2):
	"""Function that checks for the intersection of two line segments given four coordinates"""
	p_den = (start1[0] - end1[0])*(start2[1] - end2[1]) - (start1[1] - end1[1])*(start2[0] - end2[0])
	p_x = ((start1[0]*end1[1] - start1[1]*end1[0])*(start2[0] - end2[0]) - (start1[0] - end1[0])*(start2[0]*end2[1] - start2[1]*end2[0]))/p_den
	p_y = ((start1[0]*end1[1] - start1[1]*end1[0])*(start2[1] - end2[1]) - (start1[1] - end1[1])*(start2[0]*end2[1] - start2[1]*end2[0]))/p_den

	intersect = ((min([start1[0], end1[0]]) <= p_x <= max([start1[0], end1[0]])) and
                 (min([start2[0], end2[0]]) <= p_x <= max([start2[0], end2[0]])) and
                 (min([start1[1], end1[1]]) <= p_y <= max([start1[1], end1[1]])) and
                 (min([start2[1], end2[1]]) <= p_y <= max([start2[1], end2[1]])))

	return intersect, np.array([p_x, p_y])


def createRandomAngles(no_wires, kappa):
	"""Function to generate random angle in range [-pi,pi] with a Von Mises distribution"""
	dist = np.random.vonmises(0, kappa, no_wires)
	s_param_calc = np.average(2*(np.cos(dist))**2 - 1)
	return dist, s_param_calc


def two_point_resistance(val, vec, node1, node2):
	"""Calculates the resistance between 2 points using matrix solution"""
	R12 = 0
	totalnodes = len(val)
	checkvector = (1.0/np.sqrt(totalnodes))*np.ones(totalnodes)       # The trivial vector solution
	transvec = vec.T
	for i in xrange(len(val)):
		if np.allclose(checkvector, abs(transvec[i]), rtol=1e-3, atol=1e-4):
			continue        # Skips the trivial solution if found
		if val[i] <= 0:
			continue        # Skips 0 or negative eigenvalues (infinite resistance)
		R12 += (1/abs(val[i]))*(transvec[i][node1] - transvec[i][node2])**2
	return R12


def findnode(list_of_nodes, x, y):
	"""Find the closest node to a particular (x,y) coordinate"""
	xy_coord = np.array([x, y])
	distances = []
	for n in list_of_nodes:
		node_coord = np.array([n['xint'], n['yint']])
		distances.append(np.linalg.norm(xy_coord - node_coord))
	return np.argmin(distances)


# Class definition of nanowire network array object
class WireNet:
	"""Network of nanowires"""
	def __init__(self, n, l, sdl, d, sk, wres, ires, rsh, debug=False, xsort=True):
		self.n = n                     # Number of wires
		self.lav = l                   # Average length of wire
		self.lstd = sdl                # Standard deviation of wire lengths

		self.sample_dimension = d      # Size of sample

		self.angleskew = sk            # Skew of angular distribution of wires // to x-direction (kappa in von Mises distribution)

		self.wire_res = wres                # Resistance per length of wire
		self.intersect_res = ires           # Resistance of wire-to-wire interconnection
		self.sheet_res = rsh                # Sheet resistance of matrix

		self.index = np.array(range(n))     # Wire index
		self._count = -1

		self.wire_density = self.n/(self.sample_dimension**2)

		if debug:
			np.random.seed(293423)    # Debug only to get same array each time
		self._xsortbool = xsort       # If True sorts by x, if False sorts by index

		self.wirelengths = abs(np.random.normal(self.lav, self.lstd, self.n))   # Normal distribution of wire lengths
		self.startcoords = self.sample_dimension*np.random.rand(self.n, 2)      # Random start coordinates
		self.wireangles, self.s = createRandomAngles(self.n, self.angleskew)    # Von Mises distribution of wire angles
		self.endcoords = self.startcoords + np.column_stack([self.wirelengths*np.sin(self.wireangles), self.wirelengths*np.cos(self.wireangles)])

		self.allxcoords = np.column_stack([self.startcoords.T[0], self.endcoords.T[0]])    # All x-coordinates in 2xn array

	def __iter__(self):
		return self

	def next(self):
		self._count += 1
		if self._count == self.n:
			self._count = -1
			raise StopIteration
		if self._xsortbool:
			return self.sort_by_x()[self._count]
		else:
			return self.sort_by_index()[self._count]

	def areal_coverage(self, wirediameter):
		return np.sum(self.wirelengths*wirediameter)/(self.sample_dimension)**2

	def parameters(self):
		"""Prints the wire parameters"""
		print "Number of wires: {:d}".format(self.n)
		print "Wire lengths: {:.3f} with standard deviation: {:.3f}".format(self.lav, self.lstd)
		print "Sample size: {:.2f}".format(self.sample_dimension)
		print "Wire areal density: {:.2e}".format(self.wire_density)
		print "Angular dispersion: {:.3f} (calculated S-parameter of {:.3f})".format(self.angleskew, self.s)
		print "Wire resistance per length: {:.3f}".format(self.wire_res)
		print "Wire interconnect resistance: {:.3f}".format(self.intersect_res)
		print "Matrix sheet resistance: {:.3f}".format(self.sheet_res)
		return

	def sort_by_index(self):
		"""Sorts the wires by their index"""
		output = []
		for i in self.index:
			output.append((i, self.startcoords[i][0],
                              self.startcoords[i][1],
                              self.endcoords[i][0],
                              self.endcoords[i][1],
                              self.wirelengths[i]))
		return np.array(output, dtype=[('index', int), ('x1', float), ('y1', float), ('x2', float), ('y2', float), ('length', float)])

	def sort_by_x(self):
		"""Sorts the wires by the lowest x coordinate"""
		output = []
		xs_min = np.min(self.allxcoords, axis=1)           # Finds the lower of the either startcoord or endcoord
		xs_sort = np.argsort(xs_min)
		for i in self.index:
			output.append((xs_sort[i], self.startcoords[xs_sort[i]][0],
                                       self.startcoords[xs_sort[i]][1],
                                       self.endcoords[xs_sort[i]][0],
                                       self.endcoords[xs_sort[i]][1],
                                       self.wirelengths[xs_sort[i]]))
		return np.array(output, dtype=[('index', int), ('x1', float), ('y1', float), ('x2', float), ('y2', float), ('length', float)])

	def intersections(self, noprint=True, incendpoints=True):
		"""Finds intersections/nodes of wires"""
		xs_start = np.min(self.allxcoords, axis=1)         # Lowest x coordinate of wire i.e. start of wire
		xs_end = np.max(self.allxcoords, axis=1)           # Highest x coordinate of wire i.e. end of wire

		# Combined list of all x coordinates then argsorted to give indexes
		# First half of list is start coordinates, second is end coordinates
		allxs_sort = np.argsort(np.concatenate([xs_start, xs_end]))

		searchlist = []
		intersectionlist = []

		for j in xrange(2*self.n):
			wireindex = allxs_sort[j]
			if wireindex < self.n:
				# Starts of wires
				if len(searchlist) == 0:               # Adds current wire index to searchlist if empty and moves to next wire
					searchlist.append(wireindex)
					continue
				for slwireindex in searchlist:
					# Checks for intersection between wire allxs_sort[j] and all wires in searchlist
					if not noprint:
						print "Checking wires: {:3d} and {:3d}".format(wireindex, slwireindex),
					isinter, interpos = intersectCheck(self.startcoords[wireindex], self.endcoords[wireindex], self.startcoords[slwireindex], self.endcoords[slwireindex])
					if isinter:
						intersectionlist.append((interpos[0], interpos[1], wireindex, slwireindex))     # Adds intersection to list
						if not noprint:
							print "  X"
					else:
						if not noprint:
							print "   "
				searchlist.append(wireindex)           # Adds current wire index to searchlist
			else:
				# Ends of wires
				searchlist.remove(wireindex-self.n)    # Removes the current wire index from the searchlist

		# Includes wire start and end points as additional nodes if required
		if incendpoints:
			for k in xrange(self.n):
				intersectionlist.append((self.startcoords[k][0], self.startcoords[k][1], k, k))
				intersectionlist.append((self.endcoords[k][0], self.endcoords[k][1], k, k))

		return np.array(intersectionlist, dtype=[('xint', float), ('yint', float), ('wireA', int), ('wireB', int)])

	def __conductance_matrix(self):
		"""Calculates the adjacency conductance matrix"""
		intersectionlist_array = self.intersections()             # List of intersections/nodes

		no_inter = len(intersectionlist_array)                    # Total number of intersections/nodes
		conductance_ij = (1 - np.identity(no_inter))*(1.0/self.sheet_res)    # Conductance matrix with zeros on diagonal and matrix sheet resistance for other elements

		for i in xrange(self.n):
			# Finds indices of intersection list for wire i
			ionline = np.concatenate([np.where(intersectionlist_array['wireA'] == i)[0], np.where(intersectionlist_array['wireB'] == i)[0]])
			print "Wire {:3d}".format(i)
			if len(ionline) == 0:    # Checks for no intersections on the wire
				continue
			p = []                   # List of (x,y) coordinates of intersections on the wire
			for j in ionline:
				p.append(np.array([intersectionlist_array['xint'][j], intersectionlist_array['yint'][j]]))
			p = np.array(p)

			p_sort = np.argsort(p.T[0])    # Now we calculate distances by sorting p by x coordinate

			for k in xrange(len(p)-1):
				intlen = np.linalg.norm(p[p_sort[k]] - p[p_sort[k+1]])  # Distance between adjacent intersections
				if intlen == 0:
					continue
				c = 1.0/(intlen*self.wire_res + self.intersect_res)       # Conductance between intersections + resistance between wires
				print "       {:.4f} between nodes: {:4d} and {:4d}".format(c, ionline[p_sort[k]], ionline[p_sort[k+1]])
				conductance_ij[ionline[p_sort[k]]][ionline[p_sort[k+1]]] = c      # Sets relevant elements in conductance matrix
				conductance_ij[ionline[p_sort[k+1]]][ionline[p_sort[k]]] = c      # Hermitian i.e. c_ij=c_ji

		print "There are {:d} nodes from {:d} wires of which {:d} are intersections".format(no_inter, self.n, no_inter-2*self.n)
		return conductance_ij, intersectionlist_array

	def solve(self, fulloutput=False):
		"""Solves the resistor network"""
		c_ij, intersectionlist_array = self.__conductance_matrix()              # Calculates conductance matrix (adjacency)
		c_i = np.sum(c_ij, axis=1)                    # Calculates conductance matrix (degree)
		lmatrix = c_i*np.identity(len(c_i)) - c_ij    # Laplacian matrix
		print "Solving..."
		val, vec = linalg.eig(lmatrix)                # Solves Laplacian matrix
		if fulloutput:
			return np.real(val), np.real(vec), c_ij, lmatrix, intersectionlist_array
		else:
			return np.real(val), np.real(vec)

	def plot(self, node1, node2, debug=False):
		"""Plots wires and intersection points with python turtle"""
		tu.setup(width=800, height=800, startx=0, starty=0)
		tu.setworldcoordinates(-self.lav, -self.lav, self.sample_dimension+self.lav, self.sample_dimension+self.lav)
		tu.speed(0)
		tu.hideturtle()
		for i in self.index:
			if debug:
				time.sleep(2)   # Debug only
			tu.penup()
			tu.goto(self.startcoords[i][0], self.startcoords[i][1])
			tu.pendown()
			tu.goto(self.endcoords[i][0], self.endcoords[i][1])
		tu.penup()
		intersect = self.intersections(noprint=True)
		tu.goto(intersect[node1][0], intersect[node1][1])
		tu.dot(10, "blue")
		tu.goto(intersect[node2][0], intersect[node2][1])
		tu.dot(10, "blue")
		for i in intersect:
			tu.goto(i[0], i[1])
			tu.dot(4, "red")
		tu.done()
		return "Plot complete"
