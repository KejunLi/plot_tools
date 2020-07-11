#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
from configuration_for_plot import config_plot
from extraction import extract_eps, extract_filename
from sort_files import files_in_dir, sort_var_and_f

################################### Input ######################################
directory = "/home/likejun/kn1"
################################################################################
dir_f=directory

data = np.genfromtxt(dir_f, dtype=float)

x = range(len(data[:, 1]))

E_min = -839.257275
E = data[:, 1]
E = E - E_min
config_plot()

plt.plot(x, E, marker="o", markersize=2., color="tab:blue")
plt.xlabel('ionic step')
plt.ylabel('E')
plt.show()
