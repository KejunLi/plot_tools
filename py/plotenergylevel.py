#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np 
import scipy as spy
import pylab as ply
import argparse
import sys

spinup1 = [-4.0268, -2.1757, -2.0078, -1.2181, -0.5988, -0.4394, -0.1065]
spindown1 = [-3.9899, -0.9415, -0.8730, -0.3789, -0.0166, 0.2033, 0.3266]
fermi = -1.8078

for yv in spinup1:
    x = [0.2,0.5]
    y = [yv-fermi, yv-fermi]

    if yv-fermi <= 0:
        plt.plot(x,y,color='red')
    else:
        plt.plot(x,y,color='green')

for yv in spindown1:
    x = [0.5,0.8]
    y = [yv-fermi, yv-fermi]

    if yv-fermi <= 0:
        plt.plot(x,y,color='violet')
    else:
        plt.plot(x,y,color='blue')


# add label to y-axis
ply.ylabel("E [Ry]")

# set x and y range
ply.xlim([0,1])
ply.ylim([-3,2])

# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())

# show plot
#plt.show()

# save plot in a eps file
plt.savefig('spindown.eps')
