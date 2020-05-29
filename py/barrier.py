#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from constants import *
from configuration_for_plot import config_plot
from extraction import extract_etot
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_cellpara, extract_aps
from general_functions import cal_dQ, get_dQ_from_scf

################################################################################
directory = "/home/likejun/work/tibn/nk331/barrier-singlet/6x6/test"
################################################################################

(set_dir_scfout, dQ) = get_dQ_from_scf(directory)

list_dir_f = sort_var_and_f(set_dir_scfout[0])[1]
ratio = np.array(sort_var_and_f(set_dir_scfout[0])[0])
list_etot = []
for i in range(len(list_dir_f)):
    list_etot.append(float(extract_etot(list_dir_f[i], "!    total energy")[0]))


x = ratio * dQ
etot_min=min(list_etot)

y = (np.asarray(list_etot)-etot_min)

E_barrier = max(y)- y[-1]
rate_eff = np.power(10.0,12) * np.exp(-E_barrier*ev2J/(kB*T_room))
time_eff = 1.0/rate_eff
rate_eff = "{:e}".format(rate_eff, ".5f")
time_eff = "{:e}".format(time_eff, ".5f")
sys.stdout.write("\rE_barrier = {} eV".format(E_barrier))
sys.stdout.write("\rrate_eff = {} s^-1".format(rate_eff))
sys.stdout.write("\rtime_eff = {} s".format(time_eff))
sys.stdout.flush()

config_plot()
plt.plot(x, y, color="tab:blue", marker="o", markersize=5)
plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel('E (eV)')
#plt.xlim([-0.1,3.5])
#plt.ylim([-0.02,0.53])
plt.show()
