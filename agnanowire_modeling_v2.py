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

__author__ = "Jeremy Smith"
__version__ = "2.0"

# Define constants (dimensions in microns)
substratesize = 200.0
# Point to test resistance as fraction of substrate size
testpoints = substratesize*np.array([[[0.250, 0.750], [0.750, 0.250], ],
                                     [[0.250, 0.250], [0.750, 0.750], ],
                                     [[0.375, 0.625], [0.625, 0.375], ],
                                     [[0.375, 0.375], [0.625, 0.625], ],
                                     [[0.323, 0.500], [0.677, 0.500], ],
                                     [[0.500, 0.323], [0.500, 0.677]]])
nwlength = 14.0                # Nanowire length
nwlength_sd = 4.0              # Standard deviation of wire lengths
nwdiameter = 0.033             # Nanowire diameter
agresistivity = 1.59e-2                                # Ag resistivity
nwresistance = 4*agresistivity/(np.pi*nwdiameter**2)   # Nanowire resistance per unit length

nwanglekappa = 0                            # Angular distribution of wires
nwdensity = 0.020                           # Nanowires per sq micron
nwn = int(nwdensity*substratesize**2)       # Number of nanowires
nwinterres = 1.0                            # Resistance between wires
matrixrsheet = 1e8                          # Sheet resistance of matrix

# Path name for location of script
path = os.path.dirname(os.path.abspath(__file__))

nets = []
results = []

if __name__ == "__main__":
	print "\n================="
	print "Ag Nanowire Model"
	print "Jeremy Smith"
	print "=================\n"

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

	print "\n=================\n"
	nets[0].parameters()

	r_list = []   # List of resistance values
	t_list = []   # List of transmission values

	for i, n in enumerate(nets):
		eigenvalues, eigenvectors, conductance_matrix, laplacian_matrix, list_of_nodes = results[i]
		for p in testpoints:
			# Find closest nodes to the test points
			node1 = findnode(list_of_nodes, p[0][0], p[0][1])
			node2 = findnode(list_of_nodes, p[1][0], p[1][1])
			# Calculate resistance between nodes
			r = two_point_resistance(eigenvalues, eigenvectors, node1, node2)
			r_list.append(r)
			print "  R: {:.4e}".format(r)

		# Calculate transmission based on (1 - areal coverage)
		t = 100*(1 - n.areal_coverage(nwdiameter))
		t_list.append(t)
		print "  Transmission: {:.3f}".format(t)

		# Output files
		#n.sort_by_x()
		#list_of_nodes
		#laplacian_matrix
		#eigenvalues
		#eigenvectors

	print "\n=================\n"

	print "Average resistance: {:.4e} +/- {:.4e}".format(np.average(r_list), np.std(r_list)/np.sqrt(len(r_list)))
	print "Average transmission: {:.3f} +/- {:.3f}".format(np.average(t_list), np.std(t_list)/np.sqrt(len(t_list)))
