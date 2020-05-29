#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
from configuration_for_plot import config_plot, view_3d
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps, extract_cellpara
from general_functions import cal_dQ

################################### Input ######################################
directory = "/home/likejun/work/tibn/tibn_oncv_c1/11x11/nonradiative"
xlabel = "x ($\AA$)"
ylabel = "y ($\AA$)"
zlabel = "\u0394Q (amu$^{1/2}$$\AA$)"
xlim = [-13,26]
ylim = [-1,25]
zlim = [0,1.8]
title = None
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
    list_cellpara = np.asarray(extract_cellpara(dir_x))
    list_atom = extract_aps(dir_x)[0]
    list_atompos = extract_aps(dir_x)[1]
    set_atom.append(list_atom)
    set_atompos.append(np.asarray(list_atompos))

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
    set_dQ2 = []
    for j in range(len(set_atompos[i])):
        (x0, y0, z0) = np.matmul(set_atompos[i][j], list_cellpara)
        (xi, yi, zi) = np.matmul(set_atompos[i+1][j], list_cellpara)
        set_x0.append(x0)
        set_y0.append(y0)
        set_z0.append(z0)
        set_dx.append(xi-x0)
        set_dy.append(yi-y0)
        set_dz.append(zi-z0)
        dQ = cal_dQ(xi-x0, yi-y0, zi-z0, set_atom[i][j])
        set_dQ.append(dQ)
        set_dQ2.append(dQ**2)
    sys.stdout.write("\u0394Q (amu^1/2/A) = {}".format(np.sqrt(sum(set_dQ2))))
    sys.stdout.flush()

    view_3d(set_x0, set_y0, set_dQ, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
            xlim=xlim, ylim=ylim, zlim=zlim,
            title=title, view_direction="front_view")
    view_3d(set_x0, set_y0, set_dQ, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
            xlim=xlim, ylim=ylim, zlim=zlim,
            title=title, view_direction="top_view")
#    view_3d(set_x0, set_y0, set_dQ, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
#            title=title, view_direction="left_view", plot_corlorbar=True)
    break

plt.show()
