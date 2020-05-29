#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot
from extraction import extract_eps, extract_filename
from sort_files import files_in_dir, sort_var_and_f

################################### Input ######################################
directory = "/home/likejun/OSZICAR"
################################################################################
dir_f=directory
#dir_f = files_in_dir(directory, "OSZICAR")[1][0]

data = np.genfromtxt(dir_f, dtype=float)

estep = data[:, 1]
logE = np.sign(data[:, 3])*np.log(abs(data[:, 3]))

config_plot()

plt.plot(estep, logE, marker="o", markersize=2., color="tab:blue")
plt.xlabel('eletronic step')
plt.ylabel('sign(dE)*log|dE|')
plt.show()
