#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot_energy_level, plot_energy_level
from sort_files import files_in_dir
from extraction import extract_eigenenergy
from general_functions import energy_level
from constants import *

#################################### Input #####################################
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk331_nbnd1000"
vacuum = 0.180878284 # electrostatic potential in vacuum region
vbm_pristine = 1
style = 1 # 1: no label; 2: normal label; 3: specific label
reference = 1 # 1: vacuum; 2: vbm_pristine
################################################################################
"""
dir_f = files_in_dir(directory, "scf.out")[1][0]
list_E_spinu = np.asarray(extract_eigenenergy(dir_f)[0])
list_E_spind = np.asarray(extract_eigenenergy(dir_f)[1])
list_occ_spinu = np.asarray(extract_eigenenergy(dir_f)[2])
list_occ_spind = np.asarray(extract_eigenenergy(dir_f)[3])

(spinu_E, spind_E, occ_spinu, occ_spind, bands_range, vbm, cbm) = energy_level(
list_E_spinu, list_E_spind, list_occ_spinu, list_occ_spind,
lowerlimit=4, upperlimit=8)
"""

spinu_E = np.array([-6.8898, -6.8685, -6.2655, -6.1037, -3.7397, -3.6785, 3.2090, 5.0379, 3.3719])
spind_E = np.array([-6.8730, -6.6651, -6.2377, -6.0543, 3.5971, 4.5039, 3.7234, 3.7945, 5.1172])
occ_spinu = [1, 1, 1, 1, 1, 1, 0, 0, 0]
occ_spind = [1, 1, 1, 1, 0, 0, 0, 0, 0]
bands_range = [142, 143, 144, 145, 146, 147, 148, 149, 150]
spinu_E=spinu_E-1.390
spind_E=spind_E-1.390

config_plot_energy_level()
if style == 1:
    if reference == 1:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vacuum * Ry2eV)
    elif reference == 2:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vbm_pristine)
elif style == 2:
    if reference == 1:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vacuum * Ry2eV, spinu_label=None, spind_label=None)
    elif reference == 2:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vbm_pristine, spinu_label=None, spind_label=None)
else:
    if reference == 1:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vacuum * Ry2eV, spinu_label=spinu_label, spind_label=spind_label)
    elif reference == 1:
        plot_energy_level(spinu_E, spind_E, occ_spinu, occ_spind, bands_range,
            vbm_pristine, spinu_label=spinu_label, spind_label=spind_label)


plt.ylabel("E (eV)")
plt.xlim([0,1])
#plt.ylim(-8.3, 0.9) # GW eigenvalues
#plt.ylim(-8.8,0.8) # PBE0
#plt.ylim(-6.8,0.31) # PBE
# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
