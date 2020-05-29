#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from configuration_for_plot import config_plot
from extraction import extract_eps, extract_filename
from sort_files import files_in_dir, sort_var_and_f

################################### Input ######################################
directory = "/home/likejun/work/tibn/nk331/re_tibn_oncv_c1/6x6/nonradiative/yambo/data"
################################################################################

dir_f = files_in_dir(directory, "o-y.exc_amplitude_at_2")[1]
config_plot()

data = np.genfromtxt(dir_f[0], dtype=None)
#data = np.genfromtxt(i, dtype=None)
print(data)
E = data[:, 0]
eps_im = data[:, 1]
plt.plot(E, eps_im, color="tab:blue")
plt.xlabel("E (eV)")
plt.ylabel("Amplitude")
plt.show()
