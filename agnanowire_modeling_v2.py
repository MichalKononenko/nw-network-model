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
import numpy as np
from nwnet import *
import myfunctions as mf

__author__ = "Jeremy Smith"
__version__ = "2.0"

# Define constants (dimensions in microns)
substratesize = 50.0
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
nwdensity = 0.02                                                # Nanowires per sq micron
nwnumber = int(nwdensity*substratesize**2)                      # Number of nanowires
nwinterres = 10                                                 # Resistance between wires
matrixrsheet = 1.0e8                                            # Sheet resistance of matrix

# Path name for location of script
path = os.path.dirname(os.path.abspath(__file__))

summary_list = []
summary_list_header = [["filename", "nwlength", "nwstd", "nwdiameter", "agresistivity", "nwresistance",
                       "nwanglekappa", "nwdensity", "nwn", "nwinterres", "matrixrsheet",
                       "resistance", "medresistance", "resistanceerr", "transmission", "transmissionerr"]]

if __name__ == "__main__":
	print "\n================="
	print "Ag Nanowire Model"
	print "Jeremy Smith"
	print "=================\n"


	# Create nanowire network class
	n = WireNet(nwnumber, nwlength, nwlength_sd, substratesize, nwanglekappa, nwresistance, nwinterres, matrixrsheet)
	n.solve()
	node1 = findnode(n.list_of_nodes, testpoints[0][0][0], testpoints[0][0][1])
	node2 = findnode(n.list_of_nodes, testpoints[0][1][0], testpoints[0][1][1])
	n.plot(node1, node2)

	print "DONE"
