#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from constants import Ry2eV
from configuration_for_plot import config_plot

config_plot()

def neg(e, x, vbm, fac):
    y = e - (np.array(x) + vbm) * fac
    return(y)

def pos(e, x, vbm, fac):
    y = e + (np.array(x) + vbm) * fac
    return(y)


# 3x3 CB
vbm = -2.4817 - 0.242585528 * Ry2eV#pbe
cbm = 2.0481 - 0.242585528 * Ry2eV#pbe
ep = -231.0478629837 * Ry2eV # pristine
e0 = -236.4305589969 * Ry2eV
e_p1 = -236.4376285788 * Ry2eV + (0.14656387) * Ry2eV * 2
#plt.axvline(e0-e_p1-vbm, color='k', ls='dotted')
x = np.linspace(0, cbm-vbm)
e0 = e0 + x*0
e_p1 = pos(e_p1, x, vbm, 1)
#plt.plot(x, e0-ep, label="q=0 (reproduce)")
#plt.plot(x, e_p1-ep, label="q=+1 (reproduce)")


# jdftx only
vbm_jdftx = -0.213929 * Ry2eV *2
ep_0 = -116.370786 *Ry2eV *2 #energy of perfect supercell
e0_jdftx = -119.0762383942228695 * Ry2eV *2
e_p1_jdftx = (-118.1630836770944342 + 0.09000473)*Ry2eV*2
plt.axvline(e0_jdftx-e_p1_jdftx-vbm_jdftx, color='k', ls='dotted')
e0_jdftx += x*0
e_p1_jdftx = pos(e_p1_jdftx, x, vbm_jdftx, 1)
plt.plot(x, e0_jdftx-ep_0, label="q=0 (reproduced)")
plt.plot(x, e_p1_jdftx-ep_0, label="q=+1 (reproduced)")


# tutorial
vbm_o = -0.213929 * Ry2eV *2
ep_0 = -116.370786 *Ry2eV *2 #energy of perfect supercell
e0_o= -119.076238 * Ry2eV *2
e_p1_o = (-119.018475 + 0.091225)*Ry2eV*2
plt.axvline(e0_o-e_p1_o-vbm_o, color='k', ls='dotted')
e0_o = e0_o + x*0
e_p1_o = pos(e_p1_o, x, vbm_o, 1)
plt.plot(x, e0_o-ep_0, label="q=0 (tutorial)")
plt.plot(x, e_p1_o-ep_0, label="q=+1 (tutorial)")




plt.legend()
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.xlabel("$\epsilon_{fermi}-VBM$ (eV)")
plt.ylabel("Formation Energy (eV)")
plt.show()