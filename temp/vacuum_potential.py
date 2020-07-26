#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from sort_files import files_in_dir
from constants import *

################################### Input ######################################
directory = "/home/likejun/ctl/cb/cb/redo_tyler_example/3-QE/pristine"
################################################################################
plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
dir_f = files_in_dir(directory, "avgprist.out")[1]

data = np.genfromtxt(dir_f[0], dtype=float)
z = np.array(data[:, 0])
V = np.array(data[:, 1]) * Ry2eV


plt.plot(z, V, linewidth=2, color='black', label="prist")

dir_f = files_in_dir(directory, "avg.out")[1]

data = np.genfromtxt(dir_f[0], dtype=float)
z = np.array(data[:, 0])
V = np.array(data[:, 1]) * Ry2eV


plt.plot(z, V, linewidth=2, color='tab:blue', label="CB (q=+1)")

dir_f = files_in_dir(directory, "avgcharge.out")[1]

data = np.genfromtxt(dir_f[0], dtype=float)
z = np.array(data[:, 0])
V = np.array(data[:, 1]) * Ry2eV


plt.plot(z, V, linewidth=2, color='tab:red', label="CB (q=+1)")

plt.ylabel("E (eV)")
plt.xlabel("z (\u212b)")
#plt.xlim([0,1])
# remove x-axis label
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.legend()
plt.show()

# save plot in a eps file
# plt.savefig('Ti_doped_hBN.eps')
