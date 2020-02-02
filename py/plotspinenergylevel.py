#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from configuration_for_plot import config_plot
from sort_files import files_in_dir
from extraction import extract_eigenenergy

#################################### Input #####################################
directory = "/home/likejun/work/hBN/Ti/supercell_88/cs/nonradiative/relax-gs"
# Fermi level, the highest occupied level whose wavefunction is delocalized,
# is number levels below the cbm
number = 2
################################################################################

dir_f = files_in_dir(directory, "relax.out")[1][0]
l_E_spinup = extract_eigenenergy(dir_f)[0]
l_E_spindown = extract_eigenenergy(dir_f)[1]
l_occ_spinup = extract_eigenenergy(dir_f)[2]
l_occ_spindown = extract_eigenenergy(dir_f)[3]

for i, occupation in enumerate(l_occ_spinup):
    if float(occupation) == 0:
        vbm_spinup = i - 1
        cbm_spinup = i
        break
for i, occupation in enumerate(l_occ_spindown):
    if float(occupation) == 0:
        vbm_spindown = i - 1
        cbm_spindown = i
        break

if float(l_E_spinup[vbm_spinup]) >= float(l_E_spindown[vbm_spindown]):
    vbm_E = float(l_E_spinup[vbm_spinup])
    #fermi = float(l_E_spinup[vbm_spinup-number])
    print("spin up")
    print("No.{} highest occupied level {} eV".format(vbm_spinup, vbm_E))
else:
    vbm_E = float(l_E_spindown[vbm_spindown])
    #fermi = float(l_E_spindown[vbm_spindown-number])
    print("spin down")
    print("No.{} highest occupied level {} eV".format(vbm_spindown, vbm_E))

if float(l_E_spinup[cbm_spinup]) <= float(l_E_spindown[cbm_spindown]):
    cbm_E = float(l_E_spinup[cbm_spinup])
    print("spin up")
    print("No.{} lowest unoccupied state {} eV".format(cbm_spinup, cbm_E))
else:
    cbm_E = float(l_E_spindown[cbm_spindown])
    print("spin down")
    print("No.{} lowest unoccupied state {} eV".format(cbm_spindown, cbm_E))

spinup_E = l_E_spinup[range(vbm_spinup-4, vbm_spinup+6)]
spindown_E = l_E_spindown[range(vbm_spinup-4, vbm_spinup+6)]
config_plot()
# set VBM as fermi level
if float(l_E_spinup[vbm_spinup-number]) <= float(l_E_spindown[vbm_spindown-number]):
    fermi = float(l_E_spindown[vbm_spindown-number])
else:
    fermi = float(l_E_spinup[vbm_spinup-number])
x_fermi = [0,1]
y_fermi = [0,0]
plt.plot(x_fermi, y_fermi, linestyle='dashed', linewidth=0.5, color='black')

# plot energy levels of spin up and down
for yv in spinup_E:
    yv = float(yv)
    x = [0.1,0.48]
    y = [yv-fermi, yv-fermi]
    if yv-fermi <= 0:
        plt.plot(x, y, linewidth=1, color='tab:red')
        plt.text(0.05, yv-fermi, str((format(yv-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
    else:
        plt.plot(x, y, linewidth=1, color='tab:blue')
        plt.text(0.05, yv-fermi, str((format(yv-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
for yv in spindown_E:
    yv = float(yv)
    x = [0.52,0.9]
    y = [yv-fermi, yv-fermi]
    if yv-fermi <= 0:
        plt.plot(x, y, linewidth=1, color='tab:red')
        plt.text(0.95, yv-fermi, str((format(yv-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')
    else:
        plt.plot(x, y, linewidth=1, color='tab:blue')
        plt.text(0.95, yv-fermi, str((format(yv-fermi, ".2f"))), fontsize=8,
        horizontalalignment='center', verticalalignment='center')

plt.ylabel("E-E$_{Fermi}$ (eV)")
plt.xlim([0,1])
# ply.ylim([-2.5,3])

# remove x-axis label
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
