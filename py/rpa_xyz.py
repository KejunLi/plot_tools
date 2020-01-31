#!/usr/bin/env python3
# plot dielectric function (epsilon) spectrum (EPS)
import matplotlib.pyplot as plt
import numpy as np
import scipy as spy
import os
import re
from configuration_for_plot import config_plot
from extraction import extract_eps
from sort_files import files_in_dir, sort_var_and_f

################################### Input #########################################
directory = "/home/likejun/work/hBN/Ti/supercell_66/rpa/dnk_3/new_data/rps_1"
pfn = "eps" # parts of a file name
###################################################################################
config_plot()

f = files_in_dir(directory, pfn)[0]
dir_f = files_in_dir(directory, pfn)[1]
for i, file_name in enumerate(f):
    data_set = extract_eps(dir_f[i])
    #color_list = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:pink", "tab:cyan"]
    color_list = ["tab:red", "tab:green", "tab:blue"]
    marker_list = ["o", "o", "o"]
    plt.plot(data_set[0], data_set[2], marker = marker_list[i], markersize = 3, color = color_list[i], label = file_name)
plt.legend(loc = "upper left")
plt.xlabel("E (eV)")
plt.title("TiBN_6x6_nk3_rps1")
plt.ylabel("Im(${\u03B5}$)")
plt.show()
