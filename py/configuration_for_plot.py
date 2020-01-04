#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

# This function configures python plot

def config_plot():
    SMALL_SIZE = 15
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 30

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    # figure size
    plt.figure(num=1, figsize=(10,7.5), dpi=120, facecolor='w', edgecolor='k')

    # make y-axis into scientific notation, for more options go to the following website
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
    # plt.ticklabel_format(axis='y', style='plain', scilimits=(0,0), useOffset=None, useMathText=True)
