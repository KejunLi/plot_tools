#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot_energy_level, plot_energy_level
from sort_files import files_in_dir
from extraction import extract_eigenenergy
from general_functions import energy_level
from qe_post_processing import qe_post_processing
from yambo_post_processing import yambo_post_processing
from constants import *

#################################### Input #####################################
directory = "/home/likejun/work/tibn/nk331/re_tibn_oncv_c1/6x6/nonradiative/scf-cdftup1"
vacuum = 0.180697681 # electrostatic potential in vacuum region
vbm_pristine = 1
style = 1 # 1: no label; 2: normal label; 3: specific label
reference = 1 # 1: vacuum; 2: vbm_pristine
################################################################################

dir_f = files_in_dir(directory, "scf.out")[1][0]
list_E_spinu = np.asarray(extract_eigenenergy(dir_f)[0])
list_E_spind = np.asarray(extract_eigenenergy(dir_f)[1])
list_occ_spinu = np.asarray(extract_eigenenergy(dir_f)[2])
list_occ_spind = np.asarray(extract_eigenenergy(dir_f)[3])

# TiBN
#vbm_up = 147
#vbm_dn = 145
#list_E_spinu[:vbm_up] = list_E_spinu[:vbm_up] + 0.9186
#list_E_spinu[vbm_up:] = list_E_spinu[vbm_up:] + 1.1566
#list_E_spind[:vbm_dn] = list_E_spind[:vbm_dn] + 0.7925
#list_E_spind[vbm_dn:] = list_E_spind[vbm_dn:] + 1.1303

# MoBN
#vbm_up = 148
#vbm_dn = 146
#list_E_spinu[:vbm_up] = list_E_spinu[:vbm_up] + 0.9788
#list_E_spinu[vbm_up:] = list_E_spinu[vbm_up:] + 1.2422
#list_E_spind[:vbm_dn] = list_E_spind[:vbm_dn] + 0.7325
#list_E_spind[vbm_dn:] = list_E_spind[vbm_dn:] + 1.4019

(spinu_E, spind_E, occ_spinu, occ_spind, bands_range, vbm, cbm) = energy_level(
list_E_spinu, list_E_spind, list_occ_spinu, list_occ_spind,
lowerlimit=2, upperlimit=9)

"""
qe = qe_post_processing("/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/gwbse_nk221_bands1000_NGsBlkXp_5ry/scf.out")
eigenE, occ = qe.read_eigenenergies()

yambo = yambo_post_processing("/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/gwbse_nk221_bands1000_NGsBlkXp_5ry/data/o-all_Bz.qp")
correction= yambo.read_corr(6)

# corrected energy levels relative to vacuum is
# E = E_0 + correction - vacuum * Ry2eV - E_fermi
# E_fermi is read from GW report. E_0 is already output with E_fermi shifted
# to 0 as reference.
E_fermi = -1.390
spinu_E = eigenE[0, 141:154] + correction[0, :] + E_fermi
spind_E = eigenE[qe.nk, 141:154] + correction[qe.nk, :] + E_fermi
occ_spinu = occ[0, 141:154]
occ_spind = occ[qe.nk, 141:154]
bands_range = np.array(range(142, 155))
print( eigenE[qe.nk, 141:154] + correction[qe.nk, :])
print(occ_spinu)
print(occ_spind)
print(bands_range)
"""
"""
qe = qe_post_processing("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/scf.out")
eigenE, occ = qe.read_eigenenergies()
qe.read_bandgap()
yambo = yambo_post_processing("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/data/o-all_Bz.qp")
correction= yambo.read_corr(6)

# corrected energy levels relative to vacuum is
# E = E_0 + correction - vacuum * Ry2eV - E_fermi
# E_fermi is read from GW report. E_0 is already output with E_fermi shifted
# to 0 as reference.
band1 = 142
band2 = 154
spinu_E = qe.eigenE[0, band1-1:band2] + yambo.corr[0, :]
spind_E = qe.eigenE[qe.nk, band1-1:band2] + yambo.corr[yambo.nk, :]
occ_spinu = occ[0, band1-1:band2]
occ_spind = occ[qe.nk, band1-1:band2]
bands_range = np.array(range(band1, band2+1))
print(qe.eigenE[0, band1-1:band2])
print(bands_range)
"""
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
#plt.ylim(-8, 0.9) # GW eigenvalues + GW corrections
#plt.ylim(-8.3, 0.9) # PBE eigenvalues + GW corrections
#plt.ylim(-8.8,0.8) # PBE0
plt.ylim(-6.8,0.31) # PBE
# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
