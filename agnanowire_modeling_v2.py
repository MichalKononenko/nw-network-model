#! /usr/bin/env python
"""
agnanowire_modeling_v2.py
Random resitor network model for Ag nanowire network in oxide matrix
Using data from Will Scheideler

Created by Jeremy Smith on 2015-06-24
University of California, Berkeley
j-smith@ecs.berkeley.edu
"""

import os
import sys
import multiprocessing as mpr
import numpy as np
from nwnet import *
from myfunctions import *

__author__ = "Jeremy Smith"
__version__ = "2.0"

# Define constants (dimensions in microns)
substratesize = 100.0
testpoints = substratesize*np.array([[[0.250,0.750],[0.750,0.250],],
                                     [[0.250,0.250],[0.750,0.750],],
                                     [[0.375,0.625],[0.625,0.375],],
                                     [[0.375,0.375],[0.625,0.625],],
                                     [[0.323,0.500],[0.677,0.500],],
                                     [[0.500,0.323],[0.500,0.677]]])
nwlength = 14.0
nwlength_sd = 4.0
nwdiameter = 0.033
agresistivity = 1.59e-2                                # Ag resistivity
nwresistance = 4*agresistivity/(np.pi*nwdiameter**2)   # Nanowire resistance per unit length

nwanglekappa = 50
nwdensity = 0.020                           # Nanowires per sq micron
nwn = int(nwdensity*substratesize**2)       # Number of nanowires
nwinterres = 1.0
matrixrsheet = 1e8

# Path name for location of script
path = os.path.dirname(os.path.abspath(__file__))

net = []
evals = []
evecs = []
lnodes = []


def main():
	procs = 5     # Five processes running for multiprocessing
	jobs = []
	for i in xrange(procs):
		n = WireNet(nwn, nwlength, nwlength_sd, substratesize, nwanglekappa, nwresistance, nwinterres, matrixrsheet)
		net.append(n)
		n.parameters()
		eigenvalues, eigenvectors, conductance_matrix, laplacian_matrix, list_of_nodes = n.solve(fulloutput=True)
		evals.append(eigenvalues)
		evecs.append(eigenvectors)
		lnodes.append(list_of_nodes)

	for i in xrange(5):
		for p in testpoints:
			node1 = findnode(lnodes[i], p[0][0], p[0][1])
			node2 = findnode(lnodes[i], p[1][0], p[1][1])
			print two_point_resistance(evals[i], evecs[i], node1, node2)

		print "Transmission:", 100*(1 - net[i].areal_coverage(nwdiameter))

	return

if __name__ == "__main__":
	sys.exit(main())
