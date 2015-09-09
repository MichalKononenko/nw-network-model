============================================
Module for random resistor network modelling
Jeremy Smith
University of California, Berkeley
email: j-smith@eecs.berkeley.edu
============================================


Main feature is the WireNet class that creates the network of nanowires

Example of WireNet use:

from nwnet import *
net1 = WireNet(n, lav, lstd, sample_dimension, angleskew, wire_res, intersect_res, sheet_res)

In this class the following input parameters are used:
n                :    Number of wires
lav              :    Average wire length
lstd             :    Standard deviation of wire lengths
sample_dimension :    Size of sample
angleskew        :    Skew of angular distribution of wires (kappa in von Mises distribution)
wire_res         :    Resistance per length of wire
intersect_res    :    Resistance of wire-to-wire interconnection
sheet_res        :    Sheet resistance of matrix

To view these parameters as a descriptive list use:

print net1.parameters()


To calculate the conductance matrix and solve the Laplacian for the network use:

net1.solve()

This then stores the following class variables:
net1.eigenvalues      :    Eigenvalues of Laplacian matrix
net1.eigenvectors     :    Eigenvectors of Laplacian matrix
net1.laplacian        :    Laplacian matrix
net1.cmatrix          :    Conductance matrix
net1.list_of_nodes    :    List of all nodes (from intersections() function)

solve() can be run using python’s multiprocessing module as follows:

import multiprocessing
for n in nets:
    n.start()
for n in nets:
    n.join()


Then the resistance between any two points can be calculated using the helper function:

resistance12 = two_point_resistance(eigenvalues, eigenvectors, node1, node2)

(Here node1 and node2 are the indices of the nodes from list_of_nodes)


To find the node closest to a particular (x,y) coordinate use the helper function:

node = findnode(list_of_nodes, x, y)


Other useful functions are listed below:
net1.sort_by_index()    :    Returns the wires sorted by their index number
net1.sort_by_x()        :    Returns the wires sorted by their x-coordinate
net1.intersections()    :    Returns the list of nodes without solving the conductance matrix
net1.plot(node1, node2) :    Plots the wire network, nodes, and highlights 2 chosen points (node1 and node2)
net1.areal_coverage(d)  :    Calculates nanowire areal coverage as a fraction of the sample area (d = wire diameter)

Other variables:
net1.n                  :    Number of wires
net1.lav                :    Average wire length
net1.lstd               :    Standard deviation of wire lengths
net1.sample_dimension   :    Size of sample
net1.angleskew          :    Skew of angular distribution of wires (kappa in von Mises distribution)

net1.wire_res           :    Resistance per length of wire
net1.intersect_res      :    Resistance of wire-to-wire interconnection
net1.sheet_res          :    Sheet resistance of matrix

net1.s                  :    Calculated order parameter for wires (0 = unaligned, 1 = aligned)
net1.wire_density       :    Calculated areal density of wires

net1.wirelengths        :    List of wire lengths
net1.startcoords        :    List of wire start coordinates
net1.endcoords          :    List of wire end coordinates
net1.wireangles         :    List of wire angles


