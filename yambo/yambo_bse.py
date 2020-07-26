#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from plot_tools import var_color

plt.style.use("/home/likejun/work/github/plot_tools/styles/bse-publish")

dir_work = "/home/likejun/work"
dir_ti = "mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data_larger_BSEBands"
dir_mo = "mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/data"

dir_ti_y = os.path.join(dir_work, dir_ti, "o-y.alpha_q1_haydock_bse")
dir_ti_z = os.path.join(dir_work, dir_ti, "o-z.alpha_q1_haydock_bse")

dir_mo_y = os.path.join(dir_work, dir_mo, "o-y.alpha_q1_diago_bse")
dir_mo_z = os.path.join(dir_work, dir_mo, "o-z.alpha_q1_diago_bse")

dir_hbn = "pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd500_qe6.1_yambo4.4/data-12ry"
dir_bn = os.path.join(dir_work, dir_hbn, "o-x.alpha_q1_haydock_bse")

fig, ax = plt.subplots(nrows=2, ncols=2, constrained_layout=True)

data_bn = np.genfromtxt(dir_bn, dtype=float)
E_bn = data_bn[:, 0]
eps_bn = data_bn[:, 1]
eps0_bn = data_bn[:, 3]

data_ti_y = np.genfromtxt(dir_ti_y, dtype=float)
E_ti_y = data_ti_y[:, 0]
eps_ti_y = data_ti_y[:, 1]
eps0_ti_y = data_ti_y[:, 3]
data_ti_z = np.genfromtxt(dir_ti_z, dtype=float)
E_ti_z = data_ti_z[:, 0]
eps_ti_z = data_ti_z[:, 1]
eps0_ti_z = data_ti_z[:, 3]

data_mo_y = np.genfromtxt(dir_mo_y, dtype=float)
E_mo_y = data_mo_y[:, 0]
eps_mo_y = data_mo_y[:, 1]
eps0_mo_y = data_mo_y[:, 3]
data_mo_z = np.genfromtxt(dir_mo_z, dtype=float)
E_mo_z = data_mo_z[:, 0]
eps_mo_z = data_mo_z[:, 1]
eps0_mo_z = data_mo_z[:, 3]

ax[0, 0].plot(
    E_bn, eps_bn, linestyle="dashed", color=var_color("black", 0.67), 
    label="pristine hBN"
    )
ax[0, 0].plot(E_ti_y, eps_ti_y, color="tab:blue", label="along y-axis")
ax[0, 0].plot(E_ti_z, eps_ti_z, color="tab:red", label="along z-axis")
ax[0, 0].set_xlim(0.3, 2.5)
ax[0, 0].set_ylim(0, 0.3875)
#ax[0, 0].set_xlabel("E (eV)")
ax[0, 0].set_ylabel("Im(${\u03B5}$)")
ax[0, 0].legend()
ax[0, 0].xaxis.set_major_locator(plt.NullLocator())
ax[0, 0].yaxis.set_major_locator(plt.NullLocator())
coord_arrows = [(0.556, 0.125)]#, (1.427, 0.10), (1.663, 0.07), (1.986, 0.35)]
coord_trans = [(0.556, 0.18)]#, (1.3, 0.16), (1.73, 0.13), (2.25, 0.3)]
trans = [
    "1a'\u2191 \u2192 2a'\u2191", 
    #"1a''\u2191 \u2192 2a''\u2191 \n" + "...",
    #"1a''\u2191 \u2192 2a''\u2191 \n" + "1a'\u2191 \u2192 2a'\u2191 \n" + "...", 
    #"1a''\u2191 \u2192 2a''\u2191 \n" + "1a'\u2191 \u2192 2a'\u2191 \n" + "..."
    ]
"""for i in range(len(coord_arrows)):
    ax[0, 0].text(
                coord_arrows[i][0], coord_arrows[i][1], "\u2193",
                fontsize=16,
                horizontalalignment='center', 
                verticalalignment='center'
                )
for i in range(len(coord_arrows)):
    ax[0, 0].text(
                coord_trans[i][0], coord_trans[i][1], trans[i],
                fontsize=16,
                horizontalalignment='center', 
                verticalalignment='center'
                )"""


ax[0, 1].plot(
    E_bn, eps_bn, linestyle="dashed", color=var_color("black", 0.67), 
    label="pristine hBN"
    )
ax[0, 1].plot(E_ti_y, eps_ti_y, color="tab:blue", label="along y-axis")
ax[0, 1].plot(E_ti_z, eps_ti_z, color="tab:red", label="along z-axis")
ax[0, 1].set_xlim(2.5, 7.5)
ax[0, 1].set_ylim(0, 15.5)
#ax[0, 1].set_xlabel("E (eV)")
#ax[0, 1].set_ylabel("Im(${\u03B5}$)")
ax[0, 1].xaxis.set_major_locator(plt.NullLocator())
ax[0, 1].yaxis.set_major_locator(plt.NullLocator())
ax[0, 1].text(3.5, 8, "BSEBands 120-208", fontsize=20)



ax[1, 0].plot(
    E_bn, eps_bn, linestyle="dashed", color=var_color("black", 0.67), 
    label="pristine hBN"
    )
ax[1, 0].plot(E_mo_y, eps_mo_y, color="tab:blue", label="along y-axis (Mo-hBN)")
ax[1, 0].plot(E_mo_z, eps_mo_z, color="tab:red", label="along z-axis (Mo-hBN)")
ax[1, 0].set_xlim(0.3, 2.5)
ax[1, 0].set_ylim(0, 0.3875)
ax[1, 0].yaxis.set_major_locator(plt.NullLocator())
ax[1, 0].set_xlabel("E (eV)")
ax[1, 0].set_ylabel("Im(${\u03B5}$)")
#ax[1, 0].legend(loc="upper left")
coord_arrows = [(1.087, 0.08)]#, (1.100, 0.05), (1.562, 0.07)]
coord_trans = [(1.087, 0.140)]#, (1.25, 0.09), (1.73, 0.12)]
trans = [
    "2a'\u2191 \u2192 3a'\u2191", 
    #"1a'\u2191 \u2192 3a'\u2191 \n",
    #"1a''\u2191 \u2192 2a''\u2191 \n"
    ]
"""for i in range(len(coord_arrows)):
    ax[1, 0].text(
                coord_arrows[i][0], coord_arrows[i][1], "\u2193",
                fontsize=16,
                horizontalalignment='center', 
                verticalalignment='center'
                )
for i in range(len(coord_arrows)):
    ax[1, 0].text(
                coord_trans[i][0], coord_trans[i][1], trans[i],
                fontsize=16,
                horizontalalignment='center', 
                verticalalignment='center'
                )"""

ax[1, 1].plot(
    E_bn, eps_bn, linestyle="dashed", color=var_color("black", 0.67), 
    label="pristine hBN"
    )
ax[1, 1].plot(E_mo_y, eps_mo_y, color="tab:blue", label="along y-axis")
ax[1, 1].plot(E_mo_z, eps_mo_z, color="tab:red", label="along z-axis")
ax[1, 1].set_xlim(2.5, 7.5)
ax[1, 1].set_ylim(0, 15.5)
ax[1, 1].yaxis.set_major_locator(plt.NullLocator())
ax[1, 1].set_xlabel("E (eV)")
#ax[1, 1].set_ylabel("Im(${\u03B5}$)")
ax[1, 1].text(3.5, 8, "BSEBands 128-200", fontsize=20)



plt.show()
