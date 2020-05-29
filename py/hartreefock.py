#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot
from sort_files import files_in_dir
from extraction import extract_eigenenergy
from general_functions import energy_level
from constants import *
from fitting import lin_fct, best_vals_of_lin_fct


d1 = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000"
d2 = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000/pbe_scf_using_pbe0_wfc"


dir_f1 = files_in_dir(d1, "scf.out")[1][0]
E_spinu1 = np.asarray(extract_eigenenergy(dir_f1)[0])
E_spind1 = np.asarray(extract_eigenenergy(dir_f1)[1])
occ_spinu1 = np.asarray(extract_eigenenergy(dir_f1)[2])
occ_spind1 = np.asarray(extract_eigenenergy(dir_f1)[3])

dir_f2 = files_in_dir(d2, "nscf.out")[1][0]
E_spinu2 = np.asarray(extract_eigenenergy(dir_f2)[0])
E_spind2 = np.asarray(extract_eigenenergy(dir_f2)[1])
occ_spinu2 = np.asarray(extract_eigenenergy(dir_f2)[2])
occ_spind2 = np.asarray(extract_eigenenergy(dir_f2)[3])

x = np.linspace(1, len(E_spind1)+1, len(E_spind1))
HF = E_spind1-E_spind2
config_plot()

plt.plot(x, HF, "black")

best_vals = best_vals_of_lin_fct(x[145:], HF[145:])
y1 = best_vals[0] + best_vals[1]*x[145:]
sys.stdout.write("\n"+str(best_vals[0])+"\n")
sys.stdout.write(str(best_vals[1])+"\n")
plt.plot(x[145:], y1)

best_vals = best_vals_of_lin_fct(x[129:140], HF[129:140])
y2 = best_vals[0] + best_vals[1]*x[129:145]
sys.stdout.write("\n"+str(best_vals[0])+"\n")
sys.stdout.write(str(best_vals[1]))
sys.stdout.flush()
plt.plot(x[129:145], y2)

for i in range(len(x)):
    if i <= 144:
        plt.plot(x[i], HF[i], marker="o", markersize=6,
            markerfacecolor="black", color="black")
    else:
        plt.plot(x[i], HF[i], marker="o", markersize=6,
            markerfacecolor="w", color="black")


plt.ylabel("$E^{PBE0(α)}_{xc}$ - $E^{PBE}_{xc}$ (eV)")
plt.xlim([130,180])
#plt.ylim(-8.3, 0.9) # GW eigenvalues
#plt.ylim(-8.8,0.8) # PBE0
#plt.ylim(-6.8,0.31) # PBE
# remove x-axis label
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()
