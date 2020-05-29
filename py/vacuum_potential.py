#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot
from sort_files import files_in_dir
from constants import *

################################### Input ######################################
directory = "/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/scf-gs/job_vacuum"
################################################################################
config_plot()
dir_f = files_in_dir(directory, "avg.out")[1]

data = np.genfromtxt(dir_f[0], dtype=float)
z = np.array(data[:, 0])
V = np.array(data[:, 1]) * Ry2eV


plt.plot(z, V, linewidth=1, color='black')


plt.ylabel("E (eV)")
#plt.xlim([0,1])
# remove x-axis label
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
