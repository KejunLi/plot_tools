#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
sys.path.insert(0, "/home/likejun/work/github/yambo_post_processing")
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from constants import Ry2eV
import read_qe
from read_gw_output import read_gw_out



plt.style.use("/home/likejun/work/github/plot_tools/styles/sg-diag")

def plot_sg_diag(*args):
    """
    ============================================================================
    +   Default format of input:
    +   plot_sg_diagram(*args)
    +   args = [arg1, arg2, ...]
    +   arg[i] = {
            E_spinu, occ_spinu, E_spind, occ_spind, vac, 
            show_values, vbm, cbm
            }
    +
    +   arg[0]: E_spinu
    +   arg[1]: occ_spinu 
    +   arg[2]: E_spind
    +   arg[3]: occ_spind
    +   arg[4]: vac (electrostatic potential)
    +   arg[5]: show_values (0: no; 1: yes)
    +   arg[6]: vbm of pristine hBN
    +   arg[7]: cbm of pristine hBN
    +   arg[8]: plot range
    +   arg[9]: title
    +   arg[10]: number
    ============================================================================
    """
 
    fig, ax = plt.subplots(nrows=1, ncols=len(args))#, constrained_layout=True)

    for i, arg in enumerate(args):
        for j, yval in enumerate(arg[0]):
            x = [0.1, 0.45]
            y = [yval-arg[4], yval-arg[4]]
            ax[i].plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[1][j] == 1:
                    ax[i].text(
                        0.1,
                        yval - arg[4],
                        str(arg[10][j]) + "\u2b06|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    ax[i].text(
                        0.1,
                        yval - arg[4],
                        str(arg[10][j]) + "\u21e7|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[1][j] == 1:
                    ax[i].text(
                        0.25 + j%2*0.05,
                        yval-arg[4],
                        "\u2b06",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    ax[i].text(
                        0.25 + j%2*0.05,
                        yval - arg[4],
                        "\u21e7",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        for j, yval in enumerate(arg[2]):
            x = [0.55, 0.9]
            y = [yval-arg[4], yval-arg[4]]
            ax[i].plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[3][j] == 1:
                    ax[i].text(
                        0.9,
                        yval-arg[4],
                        str(arg[10][j]) + "\u2b07|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    ax[i].text(
                        0.9,
                        yval-arg[4],
                        str(arg[10][j]) + "\u21e9|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[3][j] == 1:
                    ax[i].text(
                        0.7+j%2*0.05,
                        yval-arg[4],
                        "\u2b07",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    ax[i].text(
                        0.7+j%2*0.05,
                        yval-arg[4],
                        "\u21e9",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        ax[i].set_xlim([0, 1])
        ax[i].set_ylim(arg[8])
        ax[i].fill_between(
            np.array([0, 1]),  # x range to fill
            [arg[7], arg[7]],               # lower limit of y
            #[0, 0],                         # upper limit of y
            [10, 10],                       # upper limit of y
            facecolor="tab:red",             # The fill color
            alpha=0.32,                      # Transparency of the fill
            )
        ax[i].fill_between(
            np.array([0, 1]),  # x range to fill
            [-10, -10],                     # lower limit of y
            [arg[6], arg[6]],               # upper limit of y
            facecolor="tab:blue",               # The fill color
            alpha=0.32                      # Transparency of the fill
            )
        ax[i].xaxis.set_major_locator(plt.NullLocator())
        ax[0].set_ylabel("E (eV)")
        ax[i].set_title(arg[9])

    plt.subplots_adjust(
        left=0.08, bottom=0.02, right=0.98, top=0.9, wspace=None, hspace=None
        )
    fig.savefig("foo.pdf")


show_label = 1

# TiBN
"""vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
dir_1 = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data/o-all_Bz.qp"
qe = read_qe.qe_out(dir_f)
yambo = read_gw_out(dir_1)
yambo.read_corr()
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 142:155]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_pbe0 = -7.7552
cbm_pbe0 = -0.4545
ylim = [-8.2, 2]
title = "PBE0"
number = np.arange(143, 156)
pbe = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0, ylim, title, number]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
dir_1 = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data/o-all_Bz.qp"
qe = read_qe.qe_out(dir_f)
yambo = read_gw_out(dir_1)
yambo.read_corr()
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 142:155]+yambo.corr[0,:]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]+yambo.corr[4,:]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_pbe0 = -7.8572
cbm_pbe0 = -0.2192
ylim = [-8.2, 2]
title = "GW@PBE0(read ndb.QP)"
number = np.arange(143, 156)
pbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0, ylim, title, number]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_qe.qe_out(dir_f)
qe.read_eigenenergies()
qe.eigenE[0, :qe.up_ne] += 0.9186
qe.eigenE[0, qe.up_ne:] += 1.1566
qe.eigenE[qe.nk, :qe.dn_ne] += 0.7925
qe.eigenE[qe.nk, qe.dn_ne:] += 1.1303
e_spinu = qe.eigenE[0, 142:155]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_gwpbe0 = -7.8572
cbm_gwpbe0 = -0.2192
ylim = [-8.2, 2]
title = "GW@PBE0(scissor)"
number = np.arange(143, 156)
gwpbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_gwpbe0, cbm_gwpbe0, ylim, title, number]
"""

# MoBN
vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
dir_1 = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data/o-all_Bz.qp"
qe = read_qe.qe_out(dir_f)
yambo = read_gw_out(dir_1)
yambo.read_corr()
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 142:155]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_pbe0 = -7.7552
cbm_pbe0 = -0.4545
ylim = [-8.2, 2]
title = "PBE0"
number = np.arange(143, 156)
pbe = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0, ylim, title, number]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
dir_1 = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data/o-all_Bz.qp"
qe = read_qe.qe_out(dir_f)
yambo = read_gw_out(dir_1)
yambo.read_corr()
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 142:155]+yambo.corr[0,:]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]+yambo.corr[4,:]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_pbe0 = -7.8572
cbm_pbe0 = -0.2192
ylim = [-8.2, 2]
title = "read ndb.QP"
number = np.arange(143, 156)
pbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0, ylim, title, number]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_qe.qe_out(dir_f)
qe.read_eigenenergies()
qe.eigenE[0, :qe.up_ne] += 0.9788
qe.eigenE[0, qe.up_ne:] += 1.2422
qe.eigenE[qe.nk, :qe.dn_ne] += 0.7325
qe.eigenE[qe.nk, qe.dn_ne:] += 1.4019
e_spinu = qe.eigenE[0, 142:155]
occ_spinu = qe.occ[0, 142:155]
e_spind = qe.eigenE[qe.nk, 142:155]
occ_spind = qe.occ[qe.nk, 142:155]
vbm_gwpbe0 = -7.8572
cbm_gwpbe0 = -0.2192
ylim = [-8.2, 2]
title = "scissor"
number = np.arange(143, 156)
gwpbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_gwpbe0, cbm_gwpbe0, ylim, title, number]

matplotlib.get_cachedir()

plot_sg_diag(pbe, pbe0, gwpbe0)
plt.subplots_adjust(left=0.08, bottom=0.02, right=0.98, top=0.9, wspace=None, hspace=None)
plt.show()
