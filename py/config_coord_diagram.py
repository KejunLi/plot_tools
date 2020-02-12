#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from configuration_for_plot import config_plot
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_etot, extract_aps, extract_cellpara
from fitting import quadratic_fct, best_vals_of_quadratic_fct
from general_functions import cal_dQ

################################### Input ######################################
directory = "/home/likejun/work/hBN/Ti/7x7/c1/nonradiative"
min_x = -5; max_x = 5
min_y = -0.05; max_y = 0.8
label = ["TiBN (gs)", "TiBN (ex)"] # label the two curves
label1_pos = [1,0.2]
label2_pos = [0,0.5]
title = "TiBN"
lieft_shift_E_zpl = 3.5 # shift label of ZPL
right_shift_E_rel = 1.2 # shift label of E_rel
left_arrow = 2 # move left arrow to left
right_arrow = 2 # move right arrow to right
elong_arrow = 1 # elongate the two arrow
################################################################################

# this part looks for all the scf.out files and save in the list
# for ground state and excited state, respectively.
list_dir_lin = []
list_dir_ratio = []
set_dir_scfin = []
set_dir_scfout = []
list_dir_lin = files_in_dir(directory, "lin")[1]
for dir_lin in list_dir_lin:
    list_dir_ratio.append(files_in_dir(dir_lin, "ratio-")[1])
# print(dir_ratio)
for dir_ratio in list_dir_ratio:
    list_dir_scfout_temp = []
    list_dir_scfin_temp = []
    for dir_ratio_i in dir_ratio:
        list_dir_scfin_temp.append(files_in_dir(dir_ratio_i, "scf.in")[1][0])
        list_dir_scfout_temp.append(files_in_dir(dir_ratio_i, "scf.out")[1][0])
    set_dir_scfin.append(list_dir_scfin_temp)
    set_dir_scfout.append(list_dir_scfout_temp)
# print(dir_f)

# calculate dQ
set_atom = []
set_atompos = []
# look for CELL_PARAMETERS from output
for dir_f in set_dir_scfout[0]:
    if "ratio-0.0000" in dir_f or "ratio-1.0000" in dir_f:
        list_cellpara = np.asarray(extract_cellpara(dir_f))
# look for atomic positions from scf.in
for dir_f in set_dir_scfin[0]:
    if "ratio-0.0000" in dir_f or "ratio-1.0000" in dir_f:
        list_atom = extract_aps(dir_f)[0]
        list_atompos = extract_aps(dir_f)[1]
        set_atom.append(list_atom)
        set_atompos.append(np.asarray(list_atompos))
set_dQ2 = []
for i in range(len(set_atompos)):
    for j in range(len(set_atompos[i])):
        (x0, y0, z0) = np.matmul(set_atompos[i][j], list_cellpara)
        (xi, yi, zi) = np.matmul(set_atompos[i+1][j], list_cellpara)
        dQi = cal_dQ(xi-x0, yi-y0, zi-z0, set_atom[i][j])
        set_dQ2.append(dQi**2)
    break
dQ = np.sqrt(sum(set_dQ2))
print(dQ)

# this part refines the scf.out files and extracts
# the ratio of linear extrapolation and corresponding total energies
set_etot = []
set_dQ = [] # nuclear coordinate
for i, list_dir_scfout in enumerate(set_dir_scfout):
    #print(d_f)
    list_etot = []
    nuc_coord = [] # nuclear coordinate
    sorted_list_ratio = sort_var_and_f(list_dir_scfout)[0]
    sorted_list_dir_scfout = sort_var_and_f(list_dir_scfout)[1]
    for j, dir_f in enumerate(sorted_list_dir_scfout):
        nuc_coord.append(sorted_list_ratio[j]*dQ)
        list_etot.append(extract_etot(dir_f, "!")[0])
    set_dQ.append(nuc_coord)
    set_etot.append(list_etot)
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
print("E_rel = {} eV".format(E_rel)) # The energy of gs in es geometry
print("E_abs = {} eV".format(E_abs)) # absorption
print("E_em = {} eV".format(E_em)) # emission

# this part plots
config_plot()
#axes = plt.gca()
#ylim = axes.get_ylim()
x_offset = 0.1
y_offset = 0.005

for i in range(len(set_etot)):
    if min(set_etot[i]) - min_etot == 0.0:
        best_vals = best_vals_of_quadratic_fct(set_dQ[i][0:5], set_etot[i][0:5])
        min_etot_fix = quadratic_fct(x_of_min_etot, best_vals[0],\
            best_vals[1], best_vals[2])
        sec_max_etot_fix = quadratic_fct(x_of_sec_min_etot, best_vals[0],\
            best_vals[1], best_vals[2])
        E_rel_fix = float(format(sec_max_etot_fix - min_etot_fix, ".5f"))

        plt.hlines(E_rel_fix, x_of_sec_min_etot, x_of_sec_min_etot+\
            x_offset*right_arrow, linestyles="dashed")
        plt.annotate('', xy=(x_of_sec_min_etot+x_offset*right_arrow,\
            -y_offset*elong_arrow),
            xytext=(x_of_sec_min_etot+x_offset*right_arrow,\
            E_rel_fix+y_offset*elong_arrow),
            arrowprops=dict(arrowstyle="<|-|>", color = "k"))
        plt.text(x_of_sec_min_etot+x_offset*right_shift_E_rel, E_rel_fix/2.0,
                "\u0394E$_{rel}$")
        plt.text(label1_pos[0], label1_pos[1], label[0])

        x = np.arange(min_x, max_x, 0.001)
        y = []
        for k in range(len(x)):
            y_temp = quadratic_fct(x[k], best_vals[0], best_vals[1],\
            best_vals[2]) - min_etot
            y.append(y_temp)
        plt.plot(x, y, linewidth=2, color="tab:blue")

for i in range(len(set_etot)):
    if min(set_etot[i]) - sec_min_etot == 0.0:
        best_vals = best_vals_of_quadratic_fct(set_dQ[i][7:14],\
            set_etot[i][7:14])
        sec_min_etot_fix = quadratic_fct(x_of_sec_min_etot, best_vals[0],\
            best_vals[1], best_vals[2])
        max_etot_fix = quadratic_fct(x_of_min_etot, best_vals[0],\
            best_vals[1], best_vals[2])
        E_zpl_fix = float(format(sec_min_etot_fix - min_etot_fix, ".5f"))

        plt.annotate('', xy=(x_of_min_etot-x_offset*left_arrow,\
            -y_offset*elong_arrow),
            xytext=(x_of_min_etot-x_offset*left_arrow,\
            E_zpl_fix+y_offset*elong_arrow),
            arrowprops=dict(arrowstyle="<|-|>", color = "k"))
        plt.text(x_of_min_etot-x_offset*lieft_shift_E_zpl, E_zpl_fix/2.0,\
            "\u0394E")
        plt.text(label2_pos[0], label2_pos[1], label[1])

        x = np.arange(min_x, max_x, 0.001)
        y = []
        for k in range(len(x)):
            y_temp = quadratic_fct(x[k], best_vals[0], best_vals[1],\
            best_vals[2]) - min_etot
            y.append(y_temp)
        plt.plot(x, y, linewidth=2, color="tab:red")
    for j in range(len(set_etot[i])):
        etot = set_etot[i][j]-min_etot
        if min(set_etot[i]) - min_etot == 0.0:
            plt.plot(set_dQ[i][j], etot, marker="o", markersize=6,
                markerfacecolor="w", color="tab:blue")
        elif min(set_etot[i]) - sec_min_etot == 0.0:
            plt.plot(set_dQ[i][j], etot, marker="o", markersize=6,
                markerfacecolor="w", color="tab:red")


plt.hlines(0.0, min_x, max_x, linestyles="dashed")
plt.hlines(E_zpl_fix, min_x, max_x, linestyles="dashed")
plt.vlines(x_of_min_etot, -10, 10, linestyles="dashed")
plt.vlines(x_of_sec_min_etot, -10, 10, linestyles="dashed")

E_abs_fix = float(format(max_etot_fix - min_etot_fix, ".5f"))
E_em_fix = float(format(sec_min_etot_fix - sec_max_etot_fix, ".5f"))
print("E_zpl_fix = {} eV".format(E_zpl_fix)) # ZPL
print("E_rel_fix = {} eV".format(E_rel_fix)) # The energy of gs in es geometry
print("E_abs_fix = {} eV".format(E_abs_fix)) # absorption
print("E_em_fix = {} eV".format(E_em_fix)) # emission

plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title(title)
plt.xlim([min_x,max_x])
plt.ylim([min_y,max_y])
plt.show()
