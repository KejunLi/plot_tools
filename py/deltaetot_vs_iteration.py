#!/usr/bin/env python3

##################################################################
# This program extracts data from a file that coleects total     #
# energy. The total energy data are collected by 'grep ! *out.   #
# Energy difference between two nearest data points of total     #
# energy is ploted. The data points of total energy are selected #
# to start from the second one.                                  #
##################################################################
import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from configuration_for_plot import config_plot


# read *.out to find lines with total energy and extract data from lines
def extract_etot(file_name):
    with open(file_name, "r") as f:
        lines = f.readlines()
        saved_data = []
        # unit conversion
        ry2ev = 13.6056980659
        for line in lines:
            if "!!" in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[-]\d+\.\d*", line) 
                etot = float(raw_etot[0])*ry2ev
                saved_data.append(etot)
    return(saved_data)



fname1 = "/home/KEJUNLI/work/h-BN/Mo_doped/relax_hybrid_1.out"
fname2 = "/home/KEJUNLI/work/h-BN/Mo_doped/relax_hybrid_2.out"
etot = extract_etot(fname1) + extract_etot(fname2)

# start to plot from the second data point
starting_point = 5
enume = range(starting_point -1, len(etot))
x = range(starting_point, len(etot))
y = []
for i in enume:
    if i+1 < len(etot):
        y.append(etot[i+1]-etot[i])



config_plot()
# make y-axis into scientific notation, for more options go to the following website
# https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=None, useMathText=True)

plt.plot(x, y, linewidth=1.5, marker = 'o')
ply.xlabel('Iteration')
ply.ylabel('$\Delta$E (eV)')
# change the x ticks
# plt.xticks(np.arange(min(x), max(x)+1, 2.0))

plt.show()

# save plot in a eps file
# plt.savefig('Ti_deltaetot.eps')
