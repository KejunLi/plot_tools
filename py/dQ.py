#!/usr/bin/env python3
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy as sp
import scipy.interpolate
import matplotlib.pyplot as plt
from configuration_for_plot import config_plot, view_3d
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps

################################### Input ######################################
filename = "scf.in"
directory = "/home/likejun/work/hBN/Ti/supercell_66/nonradia/my_template_yinan_structure"
# edges of unit cell
a = 15.043730685
b = 15.043730685
c = 21.167089904
#directory = "/home/likejun/work/hBN/Ti/supercell_88/nonradiative"
#a = 20.058309104
#b = 20.058309104
#c = 21.167089904
#directory = "/home/likejun/work/hBN/Ti/supercell_1010/nonradiative/cal_7"
#a = 25.072884475
#b = 25.072884475
#c = 21.167089904
view_direction = "full_view"
view_direction = "front_view"
view_direction = "top_view"
view_direction = "left_view"
################################################################################
#config_plot()
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
            set_u0 = []
            set_v0 = []
            set_w0 = []
            set_ui = []
            set_vi = []
            set_wi = []
            set_du = []
            set_dv = []
            set_dw = []
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
            set_atom_mass = []
            set_all = []
            for k in range(len(set_atompos[i][j])):
                # fractional crystal coordinates
                u0 = float(set_atompos[i][0][k][0])
                ui = float(set_atompos[i][j][k][0])
                du = ui - u0
                set_u0.append(u0)
                set_ui.append(ui)
                set_du.append(du)
                v0 = float(set_atompos[i][0][k][1])
                vi = float(set_atompos[i][j][k][1])
                dv = vi - v0
                set_v0.append(v0)
                set_vi.append(vi)
                set_dv.append(dv)
                w0 = float(set_atompos[i][0][k][2])
                wi = float(set_atompos[i][j][k][2])
                dw = wi - w0
                set_w0.append(w0)
                set_wi.append(wi)
                set_dw.append(dw)
                # convert to cartesian coordinates for hexagonal
                x0 = u0 * a + v0 * b * np.cos(np.pi*2.0/3.0)
                xi = ui * a + vi * b * np.cos(np.pi*2.0/3.0)
                dx = xi - x0
                set_x0.append(x0)
                set_xi.append(xi)
                set_dx.append(dx)
                y0 = v0 * b * np.sin(np.pi*2.0/3.0)
                yi = vi * b * np.sin(np.pi*2.0/3.0)
                dy = yi - y0
                set_y0.append(y0)
                set_yi.append(yi)
                set_dy.append(dy)
                z0 = w0 * c
                zi = wi * c
                dz = zi - z0
                set_z0.append(z0)
                set_zi.append(zi)
                set_dz.append(dz)
                if set_atom[i][j][k] == "B":
                    mass_B = 10.81
                    dQ = np.sqrt(np.power(dx,2) + np.power(dy,2) + \
                    np.power(dz,2)) * np.sqrt(mass_B)
                    set_dQ.append(dQ)
                    set_atom_mass.append(mass_B)
                elif set_atom[i][j][k] == "N":
                    mass_N = 14.01
                    dQ = np.sqrt(np.power(dx,2) + np.power(dy,2) + \
                    np.power(dz,2)) * np.sqrt(mass_N)
                    set_dQ.append(dQ)
                    set_atom_mass.append(mass_N)
                elif set_atom[i][j][k] == "Ti":
                    mass_Ti = 47.87
                    dQ = np.sqrt(np.power(dx,2) + np.power(dy,2) + \
                    np.power(dz,2)) * np.sqrt(mass_Ti)
                    set_dQ.append(dQ)
                    set_atom_mass.append(mass_Ti)
                set_all.append([x0, y0, dQ])
            if j == 0:
                view_3d(set_x0, set_y0, set_dQ, view_direction, sl_ratio[j])
            elif j == 12:
                view_3d(set_x0, set_y0, set_dQ, view_direction, sl_ratio[j])
            #fig = plt.figure(num=None, figsize=(10, 7.5), dpi=200,
            #        facecolor='w', edgecolor='k')
            #ax2 = fig.add_subplot(1,1,1)
            #ax3 = fig.add_subplot(1,1,1, projection="3d")

            # this plots dQ vs xy
            #ax3.scatter(set_x0, set_y0, set_dQ, zdir="z", s=20, c=None)
            #surf = ax3.plot_trisurf(set_x0, set_y0, set_dQ, linewidth=0.2,
            #        antialiased=True, cmap=plt.cm.viridis, alpha=0.6)

            #ax3.view_init(azim=-90, elev=0)
            #ax3.w_xaxis.line.set_lw(0.)
            #ax3.set_xticks([])
            #ax3.set_xlabel("x ($\AA$)")

            #ax3.view_init(azim=-180, elev=0)
            #ax3.w_yaxis.line.set_lw(0.)
            #ax3.set_yticks([])
            #ax3.set_ylabel("y ($\AA$)")

            #ax3.view_init(azim=-90, elev=90)
            #ax3.w_zaxis.line.set_lw(0.)
            #ax3.set_zticks([])
            #ax3.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
            #ax3.set_zlabel("z ($\AA$)")

            #ax3.set_zlim(0,0.5)
            #ax2.imshow(set_all)

            #fig.colorbar(surf)
            #fig.colorbar(surf, boundaries=np.linspace(0, 0.5))

            # this plots z vs xy
            #ax3.scatter(set_xi, set_yi, set_zi, zdir="z", s=20, c=None)
            #ax3.plot_trisurf(set_xi, set_yi, set_zi, linewidth=0.2,
            #        antialiased=True, cmap=plt.cm.Spectral)

            # this plots direction of dQ
            #ax3.scatter(set_xi, set_yi, set_zi, zdir="z", s=20, c=None)
            #ax3.quiver(set_xi, set_yi, set_zi, set_dx, set_dy, set_dz,
            #        length=0.4, linewidths=1.2, normalize=True)







            # short title for subplots
            #ax3.set_title("{}".format(sl_ratio[j]))

            # long title for single plot
            #ax3.set_title("Linear extrapolation ratio = {}".format(sl_ratio[j]))
plt.show()
