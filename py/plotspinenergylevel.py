#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np 
import scipy as spy
import pylab as ply
import argparse
import sys
import re
from configuration_for_plot import config_plot

# read and extract spin-up and spin-down energy level directly from relax*.out
def extract_data():
    if "------ SPIN UP ------------" in lines[::-1]:


'''
# read and extract spin-up and spin-down energy level data
def read_extract_spinele(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
        saved_data = []
        for line in lines:
            # Searching for positive, negative, and/or decimals, you could use [+-]?\d+(?:\.\d+)
            # add support for exponential form, try [+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?
            raw_data = re.findall(r"[+-]?\d+\.\d*", line)
            for val in raw_data:
                saved_data.append(float(val))
        return(saved_data)
'''

# directory where the file is and name of the file
def dir_and_name(directory, name):
    di = str(directory)
    na = str(name)
    return(di+na)

spinupele = dir_and_name("/home/likejun/work/Mo_doped/", "spinupele")
spindownele = dir_and_name("/home/likejun/work/Mo_doped/", "spindownele")
spinup = read_extract_spinele(spinupele)
spindown = read_extract_spinele(spindownele)


#config_plot()

# set VBM as fermi level
fermi = -4.5346
x_fermi = [0,1]
y_fermi = [0,0]
plt.plot(x_fermi, y_fermi, linestyle='dotted', linewidth=0.5, color='black')

# plot energy levels of spin up and down
for yv in spinup:
    x = [0.2,0.48]
    y = [yv-fermi, yv-fermi]
    if yv-fermi <= 0:
        plt.plot(x, y, linewidth=1, color='orangered')
    else:
        plt.plot(x, y, linewidth=1, color='forestgreen')

for yv in spindown:
    x = [0.52,0.8]
    y = [yv-fermi, yv-fermi]

    if yv-fermi <= 0:
        plt.plot(x, y, linewidth=1, color='crimson')
    else:
        plt.plot(x, y, linewidth=1, color='dodgerblue')

# add label to y-axis
ply.ylabel("E-E$_{Fermi}$ (eV)")

# set x and y range
ply.xlim([0,1])
# ply.ylim([-2.5,3])

# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())

# show plot
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
