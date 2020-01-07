#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from scipy.optimize import curve_fit
from configuration_for_plot import config_plot
from sort_files import sort_var_and_f

def extract_etot(dir_f):
    """
    read *.out to find lines with total energy and extract data from lines
    """
    with open(dir_f, "r") as f:
        lines = f.readlines()
        saved_data = []
        # unit conversion
        ry2ev = 13.6056980659
        for line in lines:
            if "!" in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line) 
                etot = float(raw_etot[0])*ry2ev
                saved_data.append(etot)
    return(saved_data)


def files_in_dir(d_dir, con_f):
    """
    collect files with satisfactory names in a designated directory
    objective files should contain con_f in their names
    """
    d_dir = d_dir + "/"
    l_f = []
    l_df = []
    for f_name in os.listdir(d_dir):
        if con_f in f_name:
            l_f.append(f_name)
            l_df.append(d_dir+f_name)
    l_f_df = [l_f, l_df]
    #print(l_f)
    #print(l_df)
    return(l_f_df)

def poly_fit(x, c1, c2, c3):
    """
    polynominal fitting
    """
    y = c1 + c2*x + c3*np.power(x,2)
    return(y)

################################### Input ############################################
directory = "/home/likejun/work/hBN/Ti/nonradia/my_template_my_structure"
#dQ = 1.136793832893554 # change of nuclear coordinate
dQ = 0.6056841293581344
min_x = -0.4
max_x = 1.6
min_y = -0.05
max_y = 0.6
label = ["TiBN (gs)", "TiBN (ex)"] # label the two curves
title = "TiBN (yinan)"
######################################################################################

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
        l_etot.append(extract_etot(s_dir_f[j])[0])
    set_dQ.append(nuc_coord)
    set_etot.append(l_etot)
    #print(min_etot)

# obtain ZPL and E_rel, and preplot
min_etot = min(min(set_etot[0]), min(set_etot[1]))
max_etot = max(max(set_etot[0]), max(set_etot[1]))
for i in range(len(set_etot)):
    if min(set_etot[i]) - min_etot > 0.0:
        sec_min_etot = min(set_etot[i])
    if max(set_etot[i]) - max_etot < 0.0:
        sec_max_etot = max(set_etot[i])
        print(sec_max_etot, max_etot)
E_zpl = float(format(sec_min_etot - min_etot, ".3f"))
E_rel = float(format(sec_max_etot - min_etot, ".3f"))
print("The ZPL is {} eV".format(E_zpl, E_rel))
print("The energy of gs in es geometry is {} eV".format(E_rel))

# this part plots
config_plot()
color_list = ["tab:blue", "tab:red"]  
axes = plt.gca()
ylim = axes.get_ylim()
x_offset = 0.1
y_offset = 0.005

# plot grids and arrows
#for i in range(len(set_etot)):
#    y_1 = min(set_etot[i])-min_etot
#    y_2 = max(set_etot[i])-min_etot
#    plt.hlines(y_1, min_x, max_x, linestyles="dashed")
#    #plt.text(max(set_dQ[0])+offset*3, y_min+offset*0.3, label[i])
#    plt.text((min_x+max_x)/1.6, y_2+x_offset*0.8*np.power(-1.0,i), label[i])
#    if min(set_etot[i]) - min_etot == 0:
#        plt.hlines(y_2, max(set_dQ[0]), max(set_dQ[0])+x_offset, linestyles="solid")

plt.hlines(0.0, min_x, max_x, linestyles="dashed")
plt.hlines(E_zpl, min_x, max_x, linestyles="dashed")
plt.vlines(min(set_dQ[0]), -10, 10, linestyles="dashed")
plt.vlines(max(set_dQ[0]), -10, 10, linestyles="dashed")
plt.hlines(E_rel, max(set_dQ[0]), max(set_dQ[0])+x_offset, linestyles="solid")

plt.annotate('', xy=(min(set_dQ[0])-x_offset*2, -y_offset), 
        xytext=(min(set_dQ[0])-x_offset*2, E_zpl+y_offset),
        arrowprops=dict(arrowstyle="<|-|>", color = "k"))
plt.annotate('', xy=(max(set_dQ[0])+x_offset, -y_offset),
        xytext=(max(set_dQ[0])+x_offset, E_rel+y_offset),
        arrowprops=dict(arrowstyle="<|-|>", color = "k"))

#plt.text((min_x+max_x)/1.6, E_rel+0.1, label[0])

plt.text(min(set_dQ[0])-x_offset*3, E_zpl/2.0, "\u0394E")
plt.text(max(set_dQ[0])+x_offset*1.2, E_rel/2.0, "\u0394E$_{rel}$")
for i in range(len(set_etot)):
    if min(set_etot[i]) - min_etot == 0:
        plt.text((min_x+max_x)/1.6, E_rel+0.1, label[i])
    elif min(set_etot[i]) - sec_min_etot == 0:
        plt.text((min_x+max_x)/1.6, E_zpl-0.1, label[i])


for i in range(len(set_etot)):
    ################### fit data #######################################################
    init_vals = [0.5, 0.5, 0.5] # for c1, c2, c3
    best_vals, covar = curve_fit(poly_fit, set_dQ[i], set_etot[i], p0=init_vals)
    print("best_vals: {}".format(best_vals))
    ####################################################################################
    x = np.arange(min_x, max_x, 0.001)
    y = []
    for k in range(len(x)):
        y_temp = poly_fit(x[k], best_vals[0], best_vals[1], best_vals[2]) - min_etot
        y.append(y_temp)
    plt.plot(x, y, linewidth=2, color=color_list[i])
    for j in range(len(set_etot[i])):
        etot = set_etot[i][j]-min_etot
        # print(etot)
        plt.plot(set_dQ[i][j], etot, marker="o", markersize=6, markerfacecolor="w", color=color_list[i])


plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title(title)
plt.xlim([min_x,max_x])
plt.ylim([min_y,max_y])
plt.show()
