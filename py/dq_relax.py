#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from configuration_for_plot import config_plot, view_3d
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps
from general_functions import cal_dQ, crystal_coord_to_cartesian_coord

################################### Input ######################################
directory = "/home/likejun/work/hBN/Ti/7x7/nonradiative"
# edges of unit cell
a = 17.5
b = 17.5
c = 15
################################################################################

# this part looks for all the scf.out files and save in the list
# for ground state and excited state, respectively.
list_dir_relax = []
list_dir_f = []
list_dir_relax = files_in_dir(directory, "relax-")[1]
for dir_relax in list_dir_relax:
    list_dir_f.append(files_in_dir(dir_relax, "relax.out")[1][0])

set_atom = []
set_atompos = []
for dir_x in list_dir_f:
    list_atom = extract_aps(dir_x)[0]
    list_atompos = extract_aps(dir_x)[1]
    set_atom.append(list_atom)
    set_atompos.append(np.asarray(list_atompos))
print(set_atompos)

for i in range(len(set_atompos)):
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
    for j in range(len(set_atompos[i])):
        u0 = float(set_atompos[i][j][0])
        ui = float(set_atompos[i+1][j][0])
        v0 = float(set_atompos[i][j][1])
        vi = float(set_atompos[i+1][j][1])
        w0 = float(set_atompos[i][j][2])
        wi = float(set_atompos[i+1][j][2])
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
        dQ = cal_dQ(dx, dy, dz, set_atom[i][j])
        set_dQ.append(dQ)

    print(DQ)
    #title = "Linear extrapolation ratio = {}".format(sl_ratio[j])
    view_3d(set_x0, set_y0, set_dQ,
            view_direction="front_view", plot_corlorbar=True)
    view_3d(set_x0, set_y0, set_dQ,
            view_direction="top_view", plot_corlorbar=True)
    view_3d(set_x0, set_y0, set_dQ,
            view_direction="left_view", plot_corlorbar=True)
    break

plt.show()
