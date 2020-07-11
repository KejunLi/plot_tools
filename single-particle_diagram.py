#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
sys.path.insert(0, "/home/likejun/work/github/yambo_post_processing")
from configuration_for_plot import config_plot_energy_level
from constants import Ry2eV
import read_output
import read_gw_output


config_plot_energy_level()

def plot_sg_diagram(*args):
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
    ============================================================================
    """
    #plt.axhline(0.0, color="k", linestyle="dashed", linewidth=0.5)
    if len(args) == 1:
        bar_width = 0.5
        pass
    else:
        num_lines = len(args) - 1
        for i in range(num_lines):
            plt.axvline(
                1.0/len(args)*(i+1), color="k", 
                linestyle="solid", linewidth=0.8
                )
    bar_width = 1.0 / len(args) / 2 

    for i, arg in enumerate(args):
        for j, yval in enumerate(arg[0]):
            x = [(i+0.1)*bar_width*2, (i+0.45)*bar_width*2]
            y = [yval-arg[4], yval-arg[4]]
            plt.plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[1][j] == 1:
                    plt.text(
                        (i+0.1)*bar_width*2,
                        yval-arg[4],
                        "\u2b06|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.1)*bar_width*2,
                        yval-arg[4],
                        "\u21e7|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[1][j] == 1:
                    plt.text(
                        (i+0.25+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u2b06",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.25+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u21e7",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        for j, yval in enumerate(arg[2]):
            x = [(i+0.55)*bar_width*2, (i+0.9)*bar_width*2]
            y = [yval-arg[4], yval-arg[4]]
            plt.plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[3][j] == 1:
                    plt.text(
                        (i+0.9)*bar_width*2,
                        yval-arg[4],
                        "\u2b07|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.9)*bar_width*2,
                        yval-arg[4],
                        "\u21e9|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[3][j] == 1:
                    plt.text(
                        (i+0.7+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u2b07",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.7+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u21e9",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        
        plt.fill_between(
            np.array([i,i+1])*bar_width*2,  # x range to fill
            [arg[7], arg[7]],               # lower limit of y
            #[0, 0],                         # upper limit of y
            [10, 10],                       # upper limit of y
            facecolor="tab:red",             # The fill color
            alpha=0.32,                      # Transparency of the fill
            )
        plt.fill_between(
            np.array([i,i+1])*bar_width*2,  # x range to fill
            [-10, -10],                     # lower limit of y
            [arg[6], arg[6]],               # upper limit of y
            facecolor="tab:blue",               # The fill color
            alpha=0.32                      # Transparency of the fill
            )


show_label = 1

# TiBN
vac = 0.180697681 * Ry2eV # electrostatic potential in vacuum region
directory = "/home/likejun/work/tibn/nk331/re_tibn_oncv_c1/6x6/nonradiative/relax-gs"
dir_f = os.path.join(directory, "relax.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 145:149]
occ_spinu = qe.occ[0, 145:149]
e_spind = qe.eigenE[qe.nk, 145:146]
occ_spind = qe.occ[qe.nk, 145:146]
vbm_pbe = -5.8262
cbm_pbe = -1.1538
pbe = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe, cbm_pbe]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 145:149]
occ_spinu = qe.occ[0, 145:149]
e_spind = qe.eigenE[qe.nk, 145:146]
occ_spind = qe.occ[qe.nk, 145:146]
vbm_pbe0 = -7.7552
cbm_pbe0 = -0.4545
pbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
qe.eigenE[0, :qe.up_ne] += 0.9186
qe.eigenE[0, qe.up_ne:] += 1.1566
qe.eigenE[qe.nk, :qe.dn_ne] += 0.7925
qe.eigenE[qe.nk, qe.dn_ne:] += 1.1303
e_spinu = qe.eigenE[0, 145:148]
occ_spinu = qe.occ[0, 145:148]
e_spind = qe.eigenE[qe.nk, 145:146]
occ_spind = qe.occ[qe.nk, 145:146]
vbm_gwpbe0 = -7.8572
cbm_gwpbe0 = -0.2192
gwpbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_gwpbe0, cbm_gwpbe0]

"""
# MoBN
vac = 0.180697681 * Ry2eV # electrostatic potential in vacuum region
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/relax-gs"
dir_f = os.path.join(directory, "relax.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 145:150]
occ_spinu = qe.occ[0, 145:150]
e_spind = qe.eigenE[qe.nk, 145:149]
occ_spind = qe.occ[qe.nk, 145:149]
vbm_pbe = -5.8262
cbm_pbe = -1.1538
pbe = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe, cbm_pbe]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
e_spinu = qe.eigenE[0, 145:150]
occ_spinu = qe.occ[0, 145:150]
e_spind = qe.eigenE[qe.nk, 145:149]
occ_spind = qe.occ[qe.nk, 145:149]
vbm_pbe0 = -7.7552
cbm_pbe0 = -0.4545
pbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_pbe0, cbm_pbe0]

vac = 0.180878284 * Ry2eV
directory = "/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4"
dir_f = os.path.join(directory, "scf.out")
qe = read_output.qe_out(dir_f)
qe.read_eigenenergies()
qe.eigenE[0, :qe.up_ne] += 0.9788
qe.eigenE[0, qe.up_ne:] += 1.2422
qe.eigenE[qe.nk, :qe.dn_ne] += 0.7325
qe.eigenE[qe.nk, qe.dn_ne:] += 1.4019
e_spinu = qe.eigenE[0, 146:150]
occ_spinu = qe.occ[0, 146:150]
e_spind = qe.eigenE[qe.nk, 145:147]
occ_spind = qe.occ[qe.nk, 145:147]
vbm_gwpbe0 = -7.8572
cbm_gwpbe0 = -0.2192
gwpbe0 = [e_spinu, occ_spinu, e_spind, occ_spind, vac, show_label, vbm_gwpbe0, cbm_gwpbe0]
"""


plot_sg_diagram(pbe, pbe0, gwpbe0)
plt.xticks([1.0/6, 1.0/2, 5.0/6], ["PBE", "PBE0", "GW@PBE0"])
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.xlim(0, 1)
plt.ylim(-8.2, 0.3)
plt.ylabel("E (eV)")
plt.show()