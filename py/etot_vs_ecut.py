#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from configuration_for_plot import config_plot

# read *.out to find lines with total energy and extract data from lines
def extract_etot(d_dir, file_name):
    destinated_file = str(d_dir) + str(file_name)
    with open(destinated_file, "r") as f:
        lines = f.readlines()
        saved_data = []
        # unit conversion
        ry2ev = 13.6056980659
        for line in lines:
            if "Final" in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line) 
                etot = float(raw_etot[0])*ry2ev
                saved_data.append(etot)
    return(saved_data)

# show ecut_*.out files in a designated directory and collect ecut and etot
def files_in_dir(d_dir):
    list_files = []
    for file_name in os.listdir(d_dir):
        if "ecut" in file_name and ".out" in file_name:
            list_files.append(file_name)
    print(list_files)
    return(list_files)

# sort out $ecut with associated Final energy
def sort_ecut_fenergy(d_dir):
    list_ecut = []
    list_etot = []
    for file_name in files_in_dir(d_dir):
        raw_ecut = re.findall(r"[+-]?\d+", file_name)
        ecut = int(raw_ecut[0])
        list_ecut.append(ecut)
        list_etot.extend(extract_etot(d_dir, file_name))
    saved_list = [list_ecut, list_etot]

    # create a sorted list of ecut and etot
    sorted_list_ecut = sorted(list_ecut)
    sorted_list_etot = []

    # sorted the list of etot according to the order of sorted_list_ecut
    for i, sorted_ecut in enumerate(sorted_list_ecut):
        for j, ecut in enumerate(list_ecut):
            if ecut == sorted_ecut:
                sorted_list_etot.append(list_etot[j])
    xy_plot = [sorted_list_ecut, sorted_list_etot]
    print(xy_plot)
    return(xy_plot)

# set the minimum of etot to be 0 as the reference energy, and rearrange the others
def rearrange_etot(d_dir):
    list_xy = sort_ecut_fenergy(d_dir)
    min_etot = min(list_xy[1])
    new_list_etot = []
    for etot in list_xy[1]:
        new_list_etot.append(etot-min_etot)
    new_list = [list_xy[0], new_list_etot]
    return(new_list)


destination = "/home/likejun/work/hBN/Ti/"
# xy = sort_ecut_fenergy(destination)
xy = rearrange_etot(destination)
x_ecut = xy[0]
y_etot = xy[1]

# fname1 = "/home/KEJUNLI/work/h-BN/Mo_doped/relax_hybrid_1.out"
# fname2 = "/home/KEJUNLI/work/h-BN/Mo_doped/relax_hybrid_2.out"
# etot = extract_etot(fname1) + extract_etot(fname2)

# start to plot from the second data point
starting_point = 1
num_points = range(starting_point-1, len(x_ecut))
x = range(starting_point, len(x_ecut)+1)
y = []
for i in num_points:
    if i <= len(x_ecut):
        y.append(y_etot[i])


config_plot()
# make y-axis into scientific notation, for more options go to the following website
# https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset=None, useMathText=True)

plt.plot(x, y, linewidth=1.5, marker = 'o')
plt.xlabel('ecut')
plt.ylabel('E$_{tot}$ (eV)')
# change the x ticks
# plt.xticks(np.arange(min(x), max(x)+1, 5.0))

plt.show()

# save plot in a eps file
# plt.savefig('Ti_etot.eps')
