#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from configuration_for_plot import config_plot
from extraction import extract_eps, extract_filename
from sort_files import files_in_dir, sort_var_and_f
from var_color import var_color


plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
"""
#directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data-scissor" # TiBN
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data" # MoBN
dir_y = os.path.join(directory, "o-y.alpha_q1_diago_bse")
dir_z = os.path.join(directory, "o-z.alpha_q1_diago_bse")

d = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd500_qe6.1_yambo4.4/data-12ry"
dir_bn = os.path.join(d, "o-x.alpha_q1_haydock_bse")

fig, ax = plt.subplots(
    num=None, figsize=(12,9), dpi=120,
    facecolor='w', edgecolor='k'
    )

data_bn = np.genfromtxt(dir_bn, dtype=float)
E_bn = data_bn[:, 0]
eps_bn = data_bn[:, 1]
eps0_bn = data_bn[:, 3]

data_y = np.genfromtxt(dir_y, dtype=float)
E_y = data_y[:, 0]
eps_y = data_y[:, 1]
eps0_y = data_y[:, 3]
data_z = np.genfromtxt(dir_z, dtype=float)
E_z = data_z[:, 0]
eps_z = data_z[:, 1]
eps0_z = data_z[:, 3]


ax.plot(E_bn, eps_bn, color=var_color("grey", 0.5), label="pristine hBN")
ax.plot(E_y, eps_y, color="tab:blue", label="along y-axis")
#ax.plot(E_y, eps0_y)
ax.plot(E_z, eps_z, color="tab:red", label="along z-axis")
#ax.plot(E_z, eps0_z)


axins = zoomed_inset_axes(ax, 4, loc="center") # zoom = 6
axins.plot(E_y, eps_y, color="tab:blue")
axins.plot(E_z, eps_z, color="tab:red")
#axins.plot(E_z, eps0_z)
axins.tick_params(axis='both', which='major', labelsize=12)


# sub region of the original image
#x1, x2, y1, y2 = 0.43, 2.23, -0.1, 0.6 # TiBN
x1, x2, y1, y2 = 0.85, 2.53, -0.1, 0.4
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

plt.xticks(visible=True)
plt.yticks(visible=True)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

ax.set_xlabel("E (eV)", fontsize=16)
ax.set_ylabel("Im(${\u03B5}$) (aribituary unit)", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlim(0.3, 8.5)
#ax.set_ylim(-0.5, 16) #TiBN
ax.set_ylim(-0.5, 12) # MoBN
ax.legend(fontsize=16)
plt.draw()
plt.show()
"""
#===============================================================================
"""
dir_f = files_in_dir("/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd500_qe6.1_yambo4.4/test", "o-x.alpha_q1_haydock_bse")[1]
data = np.genfromtxt(dir_f[0], dtype=None)
E = data[:, 0]
eps_im = data[:, 1]
eps_im1 = data[:, 3]
plt.plot(E, eps_im1, color="gray", label="GW+RPA@PBE0 (pristine hBN, read ndb.QP)")
plt.plot(E, eps_im, color="tab:blue", label="GW+BSE@PBE0 (pristine hBN, read ndb.QP)")


d = files_in_dir("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data", "o-x.alpha_q1_diago_bse")[1]
data1 = np.genfromtxt(d[0], dtype=None)
E1 = data1[:, 0]
eps_1 = data1[:, 1]
eps_2 = data1[:, 3]
#plt.plot(E1, eps_2, color="tab:purple", label="GW+IP@PBE0 (Mo-hBN, nk331, x)")
#plt.plot(E1, eps_1, color="tab:red", label="GW+BSE@PBE0 (Mo-hBN, nk331, x)")



dir_f = files_in_dir("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data", "o-y.alpha_q1_diago_bse")[1]
data = np.genfromtxt(dir_f[0], dtype=None)
E = data[:, 0]
eps_3 = data[:, 1]
eps_4 = data[:, 3]
#plt.plot(E, eps_4, color="tab:purple", label="GW+IP@PBE0 (Mo-hBN, nk331, y)")
#plt.plot(E, eps_3, color="tab:red", label="GW+BSE@PBE0 (Mo-hBN, nk331, y)")


dir_f = files_in_dir("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data", "o-z.alpha_q1_diago_bse")[1]
data = np.genfromtxt(dir_f[0], dtype=None)
E = data[:, 0]
eps_im = data[:, 1]
eps_im1 = data[:, 3]
plt.plot(E, eps_im1, color="tab:purple", label="GW+IP@PBE0 (Mo-hBN, nk331, z)")
plt.plot(E, eps_im, color="tab:red", label="GW+BSE@PBE0 (Mo-hBN, nk331, z)")



plt.legend(loc=2)
plt.xlabel("E (eV)")
plt.ylabel(r"Im($\epsilon$)")
#plt.xlim(0, 7)
#plt.xlim(0, 8)
#plt.ylim(-0.15, 0.8)
#plt.ylim(-0.5, 10.8)
plt.show()
"""
#===============================================================================

d = files_in_dir("/home/likejun/work/nbnv_66_bse_k33", "o-bse.eps_q1_haydock_bse")[1]
data1 = np.genfromtxt(d[0], dtype=None)
E1 = data1[:, 0]
eps_im1 = data1[:, 1]
eps0_im1 = data1[:, 3]
plt.plot(E1, eps_im1, label="BSE-noshift", color="k")
plt.plot(E1-0.15, eps_im1, label="BSE-shift 0.15 eV", color="tab:red")
plt.xlim(0, 2.6)
plt.xlabel("E (eV)")
plt.legend()
#plt.ylim(-0.8, 13)
plt.ylabel(r"Im($\mathrm{\epsilon}$)")
plt.show()

