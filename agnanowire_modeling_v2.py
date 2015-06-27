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

nwanglekappa = 0
nwdensity = 0.010                           # Nanowires per sq micron
nwn = int(nwdensity*substratesize**2)       # Number of nanowires
nwinterres = 1.0
matrixrsheet = 1e8

# Path name for location of script
path = os.path.dirname(os.path.abspath(__file__))

nets = []
results = []

if __name__ == "__main__":
	procs = 4                      # 4 processes running for multiprocessing (one for each nanowire network)
	for i in xrange(procs):
		# Create nanowire network class
		n = WireNet(nwn, nwlength, nwlength_sd, substratesize, nwanglekappa, nwresistance, nwinterres, matrixrsheet, queue=mpr.Queue())
		# Append the nanowire network to list
		nets.append(n)
		# Start the solve process
		n.start()

	for n in nets:
		# Append results from queue to list
		results.append(n.queue.get())
		# Wait for all processes to finish before returning to main thread
		n.join()

	for i, n in enumerate(nets):
		eigenvalues, eigenvectors, conductance_matrix, laplacian_matrix, list_of_nodes = results[i]
		for p in testpoints:
			node1 = findnode(list_of_nodes, p[0][0], p[0][1])
			node2 = findnode(list_of_nodes, p[1][0], p[1][1])
			print two_point_resistance(eigenvalues, eigenvectors, node1, node2)

		print "Transmission:", 100*(1 - n.areal_coverage(nwdiameter))



