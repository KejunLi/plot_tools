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

def view_3d(x, y, z, **kwargs):
    fig = plt.figure(num=None, figsize=(10, 7.5), dpi=120,
            facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1, projection="3d")
    ax.scatter(x, y, z, zdir="z", s=20, c=None)
    surf = ax.plot_trisurf(x, y, z, linewidth=0.2,
            antialiased=True, cmap=plt.cm.viridis, alpha=0.6)

    if "circle" in kwargs:
        alpha = np.linspace(0.0, np.pi*2, 200)
        (r, x0, y0, z0) = kwargs.get("circle")
        x1 = r * np.cos(alpha) + x0
        y1 = r * np.sin(alpha) + y0
        z1 = 0.0
        ax.plot(x1, y1, z1)
    else:
        print("Tip: if you want to show circle with defects as center, " +
            "specify 'circle = (r, x0, y0, z0)' from input")

    print("\nPlot tips\n")
    # title
    if "title" in kwargs:
        if kwargs.get("title") == None:
            print("Empty title")
        else:
            ax.set_title("{}".format(kwargs.get("title")))
    else:
        print("Tip: if you want to add title to diagram, specify " +
            "'title = {}'".format("?") + "from input")

    # about view directions and axes labels
    if "view_direction" in kwargs:
        if kwargs.get("view_direction") == "front_view":
            ax.view_init(azim=-90, elev=0)
            ax.w_yaxis.line.set_lw(0.)
            ax.set_yticks([])
            if "xlabel" in kwargs:
                if kwargs.get("xlabel") == None:
                    print("Empty xlabel")
                else:
                    ax.set_xlabel(kwargs.get("xlabel"))
            else:
                print("Tip: if you want to show xlabel, specify " +
                    "'xlabel = ?' from input")
            if "zlabel" in kwargs:
                if kwargs.get("zlabel") == None:
                    print("Empty zlabel")
                else:
                    ax.set_zlabel(kwargs.get("zlabel"))
            else:
                print("Tip: if you want to show zlabel, specify " +
                    "'zlabel = ?' from input")
        elif kwargs.get("view_direction") == "top_view":
            ax.view_init(azim=-90, elev=90)
            ax.w_zaxis.line.set_lw(0.)
            ax.set_zticks([])
            if "xlabel" in kwargs:
                if kwargs.get("xlabel") == None:
                    print("Empty xlabel")
                else:
                    ax.set_xlabel(kwargs.get("xlabel"))
            else:
                print("Tip: if you want to show xlabel, specify " +
                    "'xlabel = ?' from input")
            if "ylabel" in kwargs:
                if kwargs.get("ylabel") == None:
                    print("Empty ylabel")
                else:
                    ax.set_ylabel(kwargs.get("ylabel"))
            else:
                print("Tip: if you want to show ylabel, specify " +
                    "'ylabel = ?' from input")
        else:
            ax.view_init(azim=-180, elev=0)
            ax.w_xaxis.line.set_lw(0.)
            ax.set_xticks([])
            if "ylabel" in kwargs:
                if kwargs.get("ylabel") == None:
                    print("Empty ylabel")
                else:
                    ax.set_ylabel(kwargs.get("ylabel"))
            else:
                print("Tip: if you want to show ylabel, specify " +
                    "'ylabel = ?' from input")
            if "zlabel" in kwargs:
                if kwargs.get("zlabel") == None:
                    print("Empty zlabel")
                else:
                    ax.set_zlabel(kwargs.get("zlabel"))
            else:
                print("Tip: if you want to show zlabel, specify " +
                    "'zlabel = ?' from input")
    else:
        print("Tip: if you want to take the top view of the 3D plot, " +
            "specify 'view_direction = top_view' from input, so as " +
            "for front view and left view.")
        if "xlabel" in kwargs:
            if kwargs.get("xlabel") == None:
                print("Empty xlabel")
            else:
                ax.set_xlabel(kwargs.get("xlabel"))
        else:
            print("Tip: if you want to show xlabel, specify " +
                "'xlabel = ?' from input")
        if "ylabel" in kwargs:
            if kwargs.get("ylabel") == None:
                print("Empty ylabel")
            else:
                ax.set_ylabel(kwargs.get("ylabel"))
        else:
            print("Tip: if you want to show ylabel, specify " +
                "'ylabel = ?' from input")
        if "zlabel" in kwargs:
            if kwargs.get("zlabel") == None:
                print("Empty zlabel")
            else:
                ax.set_zlabel(kwargs.get("zlabel"))
        else:
            print("Tip: if you want to show zlabel, specify " +
                "'zlabel = ?' from input")

    # range of axes
    if "xlim" in kwargs:
        if kwargs.get("xlim") == None:
            print("Empty xlim")
        else:
            ax.set_xlim(kwargs.get("xlim"))
    else:
        print("Tip: if you want to constrain x-axis, specify "+
            "'xlim = (min,max)' from input")
    if "ylim" in kwargs:
        if kwargs.get("ylim") == None:
            print("Empty ylim")
        else:
            ax.set_ylim(kwargs.get("ylim"))
    else:
        print("Tip: if you want to constrain y-axis, specify "+
            "'ylim = (min,max)' from input")
    if "zlim" in kwargs:
        if kwargs.get("zlim") == None:
            print("Empty zlim")
        else:
            ax.set_zlim(kwargs.get("zlim"))
    else:
        print("Tip: if you want to constrain z-axis, specify "+
            "'zlim = (min,max)' from input")

    # colorbar
    if "plot_corlorbar" in kwargs:
        if kwargs.get("plot_corlorbar"):
            if "setup_boundaries" in kwargs:
                fig.colorbar(surf, boundaries=kwargs.get("setup_boundaries"))
            else:
                print("Tip: if you want to set boundaries for colorbar, " +
                    "specify 'setup_boundaries = np.linspace(min,max)' " +
                    "from input")
                fig.colorbar(surf)
    else:
        print("Tip: if you want to plot colorbar, specify " +
            "'plot_corlorbar = True' from input")
