#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from configuration_for_plot import config_plot
from sort_files import files_in_dir
from extraction import extract_eigenenergy
from general_functions import energy_level

#################################### Input #####################################
directory = "/home/likejun/work/tibn/tibn_oncv_c1/10x10/nonradiative/relax-gs"
# Fermi level, the highest occupied level whose wavefunction is delocalized,
# is num of levels below the vbm
################################################################################

dir_f = files_in_dir(directory, "relax.out")[1][0]
list_E_spinu = np.asarray(extract_eigenenergy(dir_f)[0])
list_E_spind = np.asarray(extract_eigenenergy(dir_f)[1])
list_occ_spinu = np.asarray(extract_eigenenergy(dir_f)[2])
list_occ_spind = np.asarray(extract_eigenenergy(dir_f)[3])

(spinu_E, spind_E, occ_spinu, occ_spind, vbm, cbm) = energy_level(
list_E_spinu, list_E_spind,
list_occ_spinu, list_occ_spind, lowerlimit=4, upperlimit=6)

fermi=vbm

config_plot()
plt.hlines(0.0, 0, 1, linestyles="dashed", linewidth=0.6)
# plot energy levels of spin up and down
for i, yval in enumerate(spinu_E):
    x = [0.1,0.48]
    y = [yval-fermi, yval-fermi]
    if occ_spinu[i] == 1:
        plt.plot(x, y, linewidth=1, color='tab:red')
        plt.text(0.05, yval-fermi, str((format(yval-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
    else:
        plt.plot(x, y, linewidth=1, color='tab:blue')
        plt.text(0.05, yval-fermi, str((format(yval-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
for i, yval in enumerate(spind_E):
    x = [0.52,0.9]
    y = [yval-fermi, yval-fermi]
    if occ_spind[i] == 1:
        plt.plot(x, y, linewidth=1, color='tab:red')
        plt.text(0.95, yval-fermi, str((format(yval-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
    else:
        plt.plot(x, y, linewidth=1, color='tab:blue')
        plt.text(0.95, yval-fermi, str((format(yval-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')

plt.ylabel("E-E$_{Fermi}$ (eV)")
plt.xlim([0,1])
# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
