#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
from configuration_for_plot import config_plot, view_3d
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps, extract_cellpara
from general_functions import cal_dQ, crystal_coord_to_cartesian_coord

################################### Input ######################################
directory = "/home/likejun/work/hBN/Ti/7x7/nonradiative"
xlabel = "x ($\AA$)"
ylabel = "y ($\AA$)"
zlabel = "\u0394Q (amu$^{1/2}$$\AA$)"
xlim = None
ylim = None
zlim = None
title = None
################################################################################

# looks for all the scf.out files and save in the list
# for ground state and excited state, respectively.
list_dir_lin = []
list_dir_ratio = []
set_dir_scfin = []
set_dir_scfout = []
list_dir_lin = files_in_dir(directory, "lin")[1]
for dir_lin in list_dir_lin:
    list_dir_ratio.append(files_in_dir(dir_lin, "ratio-")[1])
for dir_ratio in list_dir_ratio:
    list_dir_scfin = []
    list_dir_scfout = []
    for dir_ratio_i in dir_ratio:
        list_dir_scfin.append(files_in_dir(dir_ratio_i, "scf.in")[1][0])
        list_dir_scfout.append(files_in_dir(dir_ratio_i, "scf.out")[1][0])
    set_dir_scfin.append(list_dir_scfin)
    set_dir_scfout.append(list_dir_scfout)

# obtain CELL_PARAMETERS
for dir_scfout in set_dir_scfout[0]:
    if "ratio-0.0000" in dir_scfout or "ratio-1.0000" in dir_scfout:
        list_cellpara = np.asarray(extract_cellpara(dir_scfout))

set_atompos = []
set_atom = []
for i, list_dir_scfin in enumerate(set_dir_scfin):
    sorted_list__ratio = sort_var_and_f(list_dir_scfin)[0]
    sorted_list__dir_f = sort_var_and_f(list_dir_scfin)[1]
    ith_set_atompos = []
    ith_set_atom = []
    for j in range(len(sorted_list__dir_f)):
        list_atom = extract_aps(sorted_list__dir_f[j])[0]
        list_atompos = extract_aps(sorted_list__dir_f[j])[1]
        ith_set_atom.append(list_atom)
        ith_set_atompos.append(list_atompos)
    set_atom.append(ith_set_atom)
    set_atompos.append(ith_set_atompos)


for i in range(len(set_atompos)):
    if "gs" in set_dir_scfin[i][0]:
        for j in range(len(set_atompos[i])):
            set_x0 = []
            set_y0 = []
            set_z0 = []
            set_xi = []
            set_yi = []
            set_zi = []
            set_dx = []
            set_dy = []
            set_dz = []
            set_dQ = []
            set_dQ2 = []
            for k in range(len(set_atompos[i][j])):
                # fractional crystal coordinates
                (x0, y0, z0) = np.matmul(set_atompos[i][0][k], list_cellpara)
                (xi, yi, zi) = np.matmul(set_atompos[i][j][k], list_cellpara)
                dx = xi - x0
                set_x0.append(x0)
                set_xi.append(xi)
                set_dx.append(dx)
                dy = yi - y0
                set_y0.append(y0)
                set_yi.append(yi)
                set_dy.append(dy)
                dz = zi - z0
                set_z0.append(z0)
                set_zi.append(zi)
                set_dz.append(dz)
                dQ = cal_dQ(dx, dy, dz, set_atom[i][j][k])
                set_dQ.append(dQ)
                set_dQ2.append(dQ**2)
            if j == 12:
                sys.stdout.write("\u0394Q (amu^1/2/A) = {}"\
                .format(np.sqrt(sum(set_dQ2))))
                sys.stdout.flush()

                view_3d(set_x0, set_y0, set_dQ,
                        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                        xlim=xlim, ylim=ylim, zlim=zlim,
                        title=title, view_direction="front_view")
                view_3d(set_x0, set_y0, set_dQ,
                        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                        xlim=xlim, ylim=ylim, zlim=zlim,
                        title=title, view_direction="top_view")
#               view_3d(set_x0, set_y0, set_dQ,
#                        xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
#                        title=title, view_direction="left_view")

plt.show()
