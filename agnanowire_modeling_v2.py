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
import time
import multiprocessing as mpr
import numpy as np
from nwnet import *
import myfunctions as mf

__author__ = "Jeremy Smith"
__version__ = "2.0"

# Define constants (dimensions in microns)
substratesize = 200.0
# Point to test resistance as fraction of substrate size
testpoints = substratesize*np.array([[[0.250, 0.750], [0.750, 0.250]],
                                     [[0.250, 0.250], [0.750, 0.750]],
                                     [[0.375, 0.625], [0.625, 0.375]],
                                     [[0.375, 0.375], [0.625, 0.625]],
                                     [[0.323, 0.500], [0.677, 0.500]],
                                     [[0.500, 0.323], [0.500, 0.677]]])
nwlength = 14.0                # Nanowire length
nwlength_sd = 4.0              # Standard deviation of wire lengths
nwdiameter = 0.033             # Nanowire diameter
agresistivity = 1.59e-2                                # Ag resistivity
nwresistance = 4*agresistivity/(np.pi*nwdiameter**2)   # Nanowire resistance per unit length

nwanglekappa = 0                                                # Angular distribution of wires
nwdensity = np.array([0.040, 0.020, 0.010, 0.005])              # Nanowires per sq micron
nwnumber = (nwdensity*substratesize**2).astype(int)             # Number of nanowires
nwinterres = 1.0                                         # Resistance between wires
matrixrsheet = 1e8                                       # Sheet resistance of matrix

# Path name for location of script
path = os.path.dirname(os.path.abspath(__file__))

summary_list = []
summary_list_header = [["filename", "nwlength", "nwstd", "nwdiameter", "agresistivity", "nwresistance"
                       "nwanglekappa", "nwdensity", "nwn", "nwinterres", "matrixrsheet",
                       "resistance", "medresistance", "resistanceerr", "transmission", "transmissionerr"]]

if __name__ == "__main__":
	print "\n================="
	print "Ag Nanowire Model"
	print "Jeremy Smith"
	print "=================\n"

	# Loop for different starting parameters
	for loop, nwn in enumerate(nwnumber):
		nets = []
		results = []
		runnumber = np.random.randint(10000)
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
		for j in nets[0].parameters():
			print j

		r_list = []   # List of resistance values
		t_list = []   # List of transmission values

		for i, n in enumerate(nets):
			filename = "output_{:05d}_{:s}-{:d}".format(runnumber, n.__class__.__name__, i)
			print filename
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
			mf.dataOutputGen(filename + "_WIRES.txt", path, n.sort_by_x())
			mf.dataOutputGen(filename + "_NODES.txt", path, list_of_nodes)
			mf.dataOutputGen(filename + "_VAL.txt", path, eigenvalues)
			#mf.dataOutputGen(filename + "_VEC.txt", path, eigenvectors)

		print "\n=================\n"

		r_list = np.array(r_list)
		t_list = np.array(t_list)
		r_median = np.median(r_list)
		# Remove outliers that are an order of magnitude larger or smaller than the median
		r_list_trim = r_list[np.where((r_list < 10*r_median)&(r_list > 0.1*r_median))]

		# Compute means and standard errors
		r_av = np.mean(r_list_trim)
		r_stderr = np.std(r_list_trim)/np.sqrt(len(r_list_trim))
		t_av = np.mean(t_list)
		t_stderr = np.std(t_list)/np.sqrt(len(t_list))

		print "Median resistance: {:.4e}".format(r_median)
		print "Mean resistance: {:.4e} +/- {:.4e}".format(r_av, r_stderr)
		print "Mean transmission: {:.3f} +/- {:.3f}".format(t_av, t_stderr)

		print "\n=================\n"

		summary_list.append(["output_{:05d}".format(runnumber), nwlength, nwlength_sd, nwdiameter,
			                 agresistivity, nwresistance, nwanglekappa, nwdensity[loop], nwn, nwinterres, matrixrsheet,
			                 r_av, r_median, r_stderr, t_av, t_stderr])

	mf.dataOutputHead("SUMMARY_{:d}.txt".format(int(time.time())), path, map(list, zip(*summary_list)), summary_list_header,
		format_d="%s\t %.3f\t %.3f\t %.4f\t %.3e\t %.4e\t %.3f\t %.4f\t %i\t %.3f\t %.3e\t %.4e\t %.4e\t %.4e\t %.3f\t %.3f\n",
		format_h="%s\t")

	print "DONE"
