#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# This function configures python plot

def config_plot():
    SMALL_SIZE = 16
    MIDDLE_SIZE = 20
    BIGGER_SIZE = 30

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MIDDLE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MIDDLE_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MIDDLE_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # figure size
    plt.figure(num=None, figsize=(10,7.5), dpi=120,
                facecolor='w', edgecolor='k')

    # make y-axis into scientific notation,
    # for more options go to the following website
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
    # plt.ticklabel_format(axis='y', style='plain',
    #                       scilimits=(0,0), useOffset=None, useMathText=True)

def view_3d(set_x0, set_y0, set_dQ, **kwargs):
    fig = plt.figure(num=None, figsize=(10, 7.5), dpi=200,
            facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1, projection="3d")
    ax.scatter(set_x0, set_y0, set_dQ, zdir="z", s=20, c=None)
    surf = ax.plot_trisurf(set_x0, set_y0, set_dQ, linewidth=0.2,
            antialiased=True, cmap=plt.cm.viridis, alpha=0.6)
    ax.set_xlim(-13,26)
    ax.set_ylim(-0.5,22)
    ax.set_zlim(0,0.8)

    if "view_direction" in kwargs:
        if kwargs.get("view_direction") == "front_view":
            ax.view_init(azim=-90, elev=0)
            ax.w_yaxis.line.set_lw(0.)
            ax.set_yticks([])
            ax.set_xlabel("x ($\AA$)")
            ax.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
            #ax3.set_zlabel("z ($\AA$)")
        elif kwargs.get("view_direction") == "top_view":
            ax.view_init(azim=-90, elev=90)
            ax.w_zaxis.line.set_lw(0.)
            ax.set_zticks([])
            ax.set_xlabel("x ($\AA$)")
            ax.set_ylabel("y ($\AA$)")
        else:
            ax.view_init(azim=-180, elev=0)
            ax.w_xaxis.line.set_lw(0.)
            ax.set_xticks([])
            ax.set_ylabel("y ($\AA$)")
            ax.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
        #ax3.set_zlabel("z ($\AA$)")
    else:
        ax.set_xlabel("x ($\AA$)")
        ax.set_ylabel("y ($\AA$)")
        ax.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
        #ax3.set_zlabel("z ($\AA$)")

    if "title" in kwargs:
        ax.set_title("{}".format(kwargs.get("title")))
    else:
        print("Reminder: to add title to diagram, add 'title = {}'."\
                .format("content of title"))

    if "plot_corlorbar" in kwargs:
        if kwargs.get("plot_corlorbar"):
            if "setup_boundaries" in kwargs:
                fig.colorbar(surf, boundaries=kwargs.get("setup_boundaries"))
            else:
                print("Reminder: to set up boundaries for colorbar, " +
                    "add 'setup_boundaries = np.linspace(min,max)'.")
                fig.colorbar(surf)
    else:
        print("Reminder: to plot colorbar, add 'plot_corlorbar = True'.")
