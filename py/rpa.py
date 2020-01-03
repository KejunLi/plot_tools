#!/usr/bin/env python3

# plot dielectric function (epsilon) spectrum (EPS)

import matplotlib.pyplot as plt
import numpy as np
import scipy as spy
import os
import re
from configuration_for_plot import config_plot


def extract_data(d_dir, file_name):
    raw_data = []
    E = []
    re_eps = []
    im_eps = []
    data_set = []
    #f = open(str(d_dir)+str(file_name), 'r')
    f = str(d_dir)+str(file_name)
    E = np.loadtxt(f, usecols=(0,))
    re_eps = np.loadtxt(f, usecols=(2,))
    im_eps = np.loadtxt(f, usecols=(1,))
    data_set = [E, re_eps, im_eps, file_name]
    #lines = f.readlines()
    #for line in lines:
    #    if "#" not in line:
    #        raw_data = re.findall(r"[+-]?\d+\.\d*|[+-]?\.\d*[eE]?[+-]?\d|[+-]?\d*\.\d*[eE]?[+-]?\d", line)
    #        E.append(float(raw_data[0]))
    #        im_eps.append(float(raw_data[1]))
    #        re_eps.append(float(raw_data[2]))
    #        data_set = [E, re_eps, im_eps, file_name]
    return(data_set)

def list_files(d_dir):
    f = os.listdir(d_dir)
    file_list = []
    for file_name in f:
        #if "eps" in file_name:
        if "eps" in file_name and "q2" in file_name:
            file_list.append(file_name)
    print(file_list)
    return(file_list)

def plot_pen(d_dir):
    f = list_files(d_dir)
    f = sorted(f)
    for i, file_name in enumerate(f):
        data_set = extract_data(d_dir, file_name)
        color_list = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:pink", "tab:cyan"]
        marker_list = ["o", "^", "", "", "", "", "", "", "", "", ""]
        plt.plot(data_set[0], data_set[2], marker = marker_list[i], markersize = 4-2*i, color = color_list[i+1], label = data_set[3])
    plt.legend(loc = "upper left")
    plt.xlabel("E (eV)")
    plt.ylabel("Im(${\u03B5}$)")
    plt.show()

config_plot()
destination = "/home/likejun/work/hBN/Mo/rpa/dnk_5/"
plot_pen(destination)
