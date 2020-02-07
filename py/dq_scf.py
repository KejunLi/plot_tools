#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from configuration_for_plot import config_plot, view_3d
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps
from general_functions import cal_dQ, crystal_coord_to_cartesian_coord

################################### Input ######################################
filename = "scf.in"
directory = "/home/likejun/work/hBN/Ti/7x7/nonradiative"
# edges of unit cell
a = 17.5
b = 17.5
c = 15
################################################################################

# this part looks for all the scf.out files and save in the list
# for ground state and excited state, respectively.
l_dir_lin = []
l_dir_ratio = []
l_dir_f = []
set_dir_f = []
l_dir_lin = files_in_dir(directory, "lin")[1]
for dir_lin in l_dir_lin:
    l_dir_ratio.append(files_in_dir(dir_lin, "ratio-")[1])
for dir_ratio in l_dir_ratio:
    l_dir_f_temp = []
    for dir_ratio_i in dir_ratio:
        l_dir_f_temp.append(files_in_dir(dir_ratio_i, filename)[1][0])
    set_dir_f.append(l_dir_f_temp)


set_atompos = []
set_atom = []
for i, l_dir_f in enumerate(set_dir_f):
    sl_ratio = sort_var_and_f(l_dir_f)[0]
    sl_dir_f = sort_var_and_f(l_dir_f)[1]
    ith_set_atompos = []
    ith_set_atom = []
    for j in range(len(sl_dir_f)):
        l_atom = extract_aps(sl_dir_f[j])[0]
        l_atompos = extract_aps(sl_dir_f[j])[1]
        ith_set_atom.append(l_atom)
        ith_set_atompos.append(l_atompos)
    set_atom.append(ith_set_atom)
    set_atompos.append(ith_set_atompos)


for i in range(len(set_atompos)):
    if "gs" in set_dir_f[i][0]:
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
            DQ = 0
            for k in range(len(set_atompos[i][j])):
                # fractional crystal coordinates
                u0 = float(set_atompos[i][0][k][0])
                ui = float(set_atompos[i][j][k][0])
                v0 = float(set_atompos[i][0][k][1])
                vi = float(set_atompos[i][j][k][1])
                w0 = float(set_atompos[i][0][k][2])
                wi = float(set_atompos[i][j][k][2])
                x0 = crystal_coord_to_cartesian_coord(a, b, c, u0, v0, w0)[0]
                y0 = crystal_coord_to_cartesian_coord(a, b, c, u0, v0, w0)[1]
                z0 = crystal_coord_to_cartesian_coord(a, b, c, u0, v0, w0)[2]
                xi = crystal_coord_to_cartesian_coord(a, b, c, ui, vi, wi)[0]
                yi = crystal_coord_to_cartesian_coord(a, b, c, ui, vi, wi)[1]
                zi = crystal_coord_to_cartesian_coord(a, b, c, ui, vi, wi)[2]
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
            if j == 12:
                print(DQ)
                title = "Linear extrapolation ratio = {}".format(sl_ratio[j])
                view_3d(set_x0, set_y0, set_dQ, title=title,
                        view_direction="front_view", plot_corlorbar=True)
                view_3d(set_x0, set_y0, set_dQ, title=title,
                        view_direction="top_view", plot_corlorbar=True)
                view_3d(set_x0, set_y0, set_dQ, title=title,
                        view_direction="left_view", plot_corlorbar=True)

plt.show()
