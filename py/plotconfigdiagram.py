#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from scipy.optimize import curve_fit
from configuration_for_plot import config_plot
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_etot
from fitting import poly_fct, best_vals_of_poly_fct

################################### Input ####################################################
directory = "/home/KEJUNLI/work/hBN/Ti/supercell_1010/nonradiative/cal_6"
dQ = 1.2911782126728621 # change of nuclear coordinate
min_x = -0.4; max_x = 2.0
min_y = -0.05; max_y = 0.6
label = ["TiBN (gs)", "TiBN (ex)"] # label the two curves
title = "TiBN_10x10 (kejun)"
##############################################################################################
######## Shift the labels of ZPL and E_rel ##############
############ usually don't need to change ###############
left_shift_pos_E_zpl = 3
right_shift_pos_E_rel = 1.2
#########################################################

# this part looks for all the scf.out files and save in the list for ground state and excited state, respectively.
dir_lin = []
dir_ratio = []
dir_f = []
dir_lin = files_in_dir(directory, "lin")[1]
for d_lin in dir_lin:
    dir_ratio.append(files_in_dir(d_lin, "ratio-")[1])
# print(dir_ratio)
for d_ratio in dir_ratio:
    dir_f_temp = []
    for d_rat in d_ratio:
        dir_f_temp.append(files_in_dir(d_rat, "scf.out")[1][0])
    dir_f.append(dir_f_temp)
# print(dir_f)

# this part refines the scf.out files and extracts the ratio of linear extrapolation and corresponding total energies
set_etot = []
set_dQ = [] # nuclear coordinate
for i, d_f in enumerate(dir_f):
    #print(d_f)
    l_etot = []
    nuc_coord = [] # nuclear coordinate
    s_ratio = sort_var_and_f(d_f)[0]
    s_dir_f = sort_var_and_f(d_f)[1]
    #print(s_ratio, s_dir_f)
    for j in range(len(s_dir_f)):
        nuc_coord.append(s_ratio[j]*dQ)
        l_etot.append(extract_etot(s_dir_f[j], "!")[0])
    set_dQ.append(nuc_coord)
    set_etot.append(l_etot)
    #print(min_etot)

# obtain ZPL, E_rel, E_abs and E_em, and prepare for plotting
min_etot = min(min(set_etot[0]), min(set_etot[1]))
max_etot = max(max(set_etot[0]), max(set_etot[1]))
# the minimum total energy in the first excited state
sec_min_etot = None 
# the maximum total energy of excited state in the ground state geometry
# where the ratio of nuclear coordinate change is 1.00
sec_max_etot = None 
# the 1D coordinate of the minimum total energy
x_of_min_etot = None
# the 1D coordinate of the second minimum total energy
x_of_sec_min_etot = None
for i in range(len(set_etot)):
    if min(set_etot[i]) - min_etot > 0.0:
        sec_min_etot = min(set_etot[i])
        for j in range(len(set_etot[i])):
            if set_etot[i][j] - sec_min_etot == 0.0:
                x_of_sec_min_etot = set_dQ[i][j]
for i in range(len(set_etot)):  
    if min(set_etot[i]) - min_etot == 0.0:
        for j in range(len(set_etot[i])):
            if set_etot[i][j] - min_etot == 0.0:
                x_of_min_etot = set_dQ[i][j]
            if set_dQ[i][j] - x_of_sec_min_etot == 0.0:
                sec_max_etot = set_etot[i][j]


E_zpl = float(format(sec_min_etot - min_etot, ".5f"))
E_rel = float(format(sec_max_etot - min_etot, ".5f"))
E_abs = float(format(max_etot - min_etot, ".5f"))
E_em = float(format(sec_min_etot - sec_max_etot, ".5f"))
print("E_zpl = {} eV".format(E_zpl)) # ZPL
print("E_rel {} eV".format(E_rel)) # The energy of gs in es geometry
print("E_abs = {} eV".format(E_abs)) # absorption
print("E_em = {} eV".format(E_em)) # emission

# this part plots
config_plot()
#axes = plt.gca()
#ylim = axes.get_ylim()
x_offset = 0.1
y_offset = 0.005

plt.hlines(0.0, min_x, max_x, linestyles="dashed")
plt.hlines(E_zpl, min_x, max_x, linestyles="dashed")
plt.vlines(x_of_min_etot, -10, 10, linestyles="dashed")
plt.vlines(x_of_sec_min_etot, -10, 10, linestyles="dashed")
plt.hlines(E_rel, x_of_sec_min_etot, x_of_sec_min_etot+x_offset, linestyles="dashed")

plt.annotate('', xy=(x_of_min_etot-x_offset*2, -y_offset), 
        xytext=(x_of_min_etot-x_offset*2, E_zpl+y_offset),
        arrowprops=dict(arrowstyle="<|-|>", color = "k"))
plt.annotate('', xy=(x_of_sec_min_etot+x_offset, -y_offset),
        xytext=(x_of_sec_min_etot+x_offset, E_rel+y_offset),
        arrowprops=dict(arrowstyle="<|-|>", color = "k"))

plt.text(x_of_min_etot-x_offset*left_shift_pos_E_zpl, E_zpl/2.0, "\u0394E")
plt.text(x_of_sec_min_etot+x_offset*right_shift_pos_E_rel, E_rel/2.0, "\u0394E$_{rel}$")
for i in range(len(set_etot)):
    if min(set_etot[i]) - min_etot == 0.0:
        plt.text((min_x+max_x)/1.6, E_rel+0.1, label[i])
    elif min(set_etot[i]) - sec_min_etot == 0.0:
        plt.text((min_x+max_x)/1.6, E_zpl-0.1, label[i])

for i in range(len(set_etot)):
    best_vals = best_vals_of_poly_fct(set_dQ[i], set_etot[i])
    x = np.arange(min_x, max_x, 0.001)
    y = []
    for k in range(len(x)):
        y_temp = poly_fct(x[k], best_vals[0], best_vals[1], best_vals[2]) - min_etot
        y.append(y_temp)
    if min(set_etot[i]) - min_etot == 0.0:
        plt.plot(x, y, linewidth=2, color="tab:blue")
    elif min(set_etot[i]) - sec_min_etot == 0.0:
        plt.plot(x, y, linewidth=2, color="tab:red")
    for j in range(len(set_etot[i])):
        etot = set_etot[i][j]-min_etot
        if min(set_etot[i]) - min_etot == 0.0:
            plt.plot(set_dQ[i][j], etot, marker="o", markersize=6, markerfacecolor="w", color="tab:blue")
        elif min(set_etot[i]) - sec_min_etot == 0.0:
            plt.plot(set_dQ[i][j], etot, marker="o", markersize=6, markerfacecolor="w", color="tab:red")

plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title(title)
plt.xlim([min_x,max_x])
plt.ylim([min_y,max_y])
plt.show()
