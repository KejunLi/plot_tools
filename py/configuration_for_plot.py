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
    plt.figure(num=None, figsize=(10,7.5), dpi=120, facecolor='w', edgecolor='k')

    # make y-axis into scientific notation, for more options go to the following website
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
    # plt.ticklabel_format(axis='y', style='plain', scilimits=(0,0), useOffset=None, useMathText=True)

def view_3d(set_x0, set_y0, set_dQ, view_direction, title):
    fig = plt.figure(num=None, figsize=(10, 7.5), dpi=200,
            facecolor='w', edgecolor='k')
    ax3 = fig.add_subplot(1,1,1, projection="3d")
    ax3.scatter(set_x0, set_y0, set_dQ, zdir="z", s=20, c=None)
    surf = ax3.plot_trisurf(set_x0, set_y0, set_dQ, linewidth=0.2,
            antialiased=True, cmap=plt.cm.viridis, alpha=0.6)
    if view_direction == "full_view":
        ax3.set_xlabel("x ($\AA$)")
        ax3.set_ylabel("y ($\AA$)")
        ax3.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
        #ax3.set_zlabel("z ($\AA$)")
    elif view_direction == "front_view":
        ax3.view_init(azim=-90, elev=0)
        ax3.w_yaxis.line.set_lw(0.)
        ax3.set_yticks([])
        ax3.set_xlabel("x ($\AA$)")
        ax3.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
        #ax3.set_zlabel("z ($\AA$)")
    elif view_direction == "top_view":
        ax3.view_init(azim=-90, elev=90)
        ax3.w_zaxis.line.set_lw(0.)
        ax3.set_zticks([])
        ax3.set_xlabel("x ($\AA$)")
        ax3.set_ylabel("y ($\AA$)")
    else:
        ax3.view_init(azim=-180, elev=0)
        ax3.w_xaxis.line.set_lw(0.)
        ax3.set_xticks([])
        ax3.set_ylabel("y ($\AA$)")
        ax3.set_zlabel("\u0394Q (amu$^{1/2}$$\AA$)")
        #ax3.set_zlabel("z ($\AA$)")

    ax3.set_xlim(-13,26)
    ax3.set_ylim(-0.5,22)
    ax3.set_zlim(0,0.8)
    #ax2.imshow(set_all)
    fig.colorbar(surf)
    ax3.set_title("Linear extrapolation ratio = {}".format(title))

    #fig.colorbar(surf, boundaries=np.linspace(0, 0.5))
