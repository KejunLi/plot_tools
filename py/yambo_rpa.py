#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from configuration_for_plot import config_plot
from extraction import extract_eps, extract_filename
from sort_files import files_in_dir, sort_var_and_f

################################### Input ######################################
directory = "/home/likejun/work/tibn/tibn_oncv_c1/6x6/nonradiative/yambo/rpa_using_yambo4.1/data"
title = "6x6"
label = "(100)"
################################################################################

dir_f = files_in_dir(directory, "o-rpa.eps")[1][0]

data = np.genfromtxt(dir_f, dtype=None)
E = data[:, 0]
eps_im = data[:, 1]

config_plot()

plt.plot(E, eps_im, marker="o", markersize=1., color="tab:blue", label=label)
plt.legend(loc = "upper left")
plt.xlabel("E (eV)")
plt.title(title)
plt.ylabel("Im(${\u03B5}$)")
plt.show()
