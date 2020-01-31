#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from sort_files import sort_var_and_f
from extraction import extract_filename, extract_eps
from var_color import var_color
from configuration_for_plot import config_plot

def list_files(d_dir):
    f = os.listdir(d_dir)
    l_f = []
    for f_name in f:
        #if "eps" in f_name:
        if "z.eps" in f_name:
            l_f.append(f_name)
    print(l_f)
    return(l_f)

def plot_pen(d_dir):
    l_var_f = list_files(d_dir)
    sl_var = sort_var_and_f(l_var_f)[0]
    sl_f = sort_var_and_f(l_var_f)[1]
    # print(sl_var, sl_f)
    sl_rf = []
    for f_name in sl_f:
        rf_name = extract_filename(f_name, "q")
        sl_rf.append(rf_name)
    for i, f_name in enumerate(sl_f):
        data_set = extract_eps(d_dir, f_name)
        # color_list = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:pink", "tab:cyan"]
        # marker_list = ["", "", "", "", "", "", "", "", "", "", ""]
        # plt.plot(data_set[0], data_set[2], marker = marker_list[i], markersize = 2, color = color_list[i+1], label = data_set[3])
        # plt.plot(data_set[0], data_set[2], marker = "", markersize = 2, color = var_color("tab:red", sl_var[i]*0.1), label = data_set[3])
        plt.plot(data_set[0], data_set[2], marker = "", markersize = 2, color = var_color("tab:blue", sl_var[i]*0.3), label = sl_rf[i])
    plt.legend(loc = "upper left")
    plt.xlabel("E (eV)")
    plt.ylabel("Im(${\u03B5}$)")
    plt.title("RPA")
    plt.show()

config_plot()
destination = "/home/likejun/work/hBN/Mo/rpa/dnk_3/"
plot_pen(destination)
