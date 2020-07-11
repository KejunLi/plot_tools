#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
from fit_functions import best_vals_of_quadratic_fct, quadratic_fct
from constants import *
from sort_files import files_in_dir


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
    #plt.figure(num=None, figsize=(16,12), dpi=120,
    #            facecolor='w', edgecolor='k')
    plt.figure(num=None, figsize=(12,9), dpi=120,
                facecolor='w', edgecolor='k')
    # make y-axis into scientific notation,
    # for more options go to the following website
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.ticklabel_format.html
    plt.ticklabel_format(axis='y', style='plain',
                           scilimits=(0,0), useOffset=False, useMathText=True)


def config_plot_energy_level():
    SMALL_SIZE = 20
    MIDDLE_SIZE = 16
    BIGGER_SIZE = 30

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MIDDLE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MIDDLE_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MIDDLE_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.figure(num=None, figsize=(12,9), dpi=120,
                facecolor='w', edgecolor='k')


def plot_sg_diagram(*args):
    """
    ============================================================================
    +   Default format of input:
    +   plot_sg_diagram(*args)
    +   args = [arg1, arg2, ...]
    +   arg[i] = {
            E_spinu, occ_spinu, E_spind, occ_spind, vac, 
            show_values, vbm, cbm
            }
    +
    +   arg[0]: E_spinu
    +   arg[1]: occ_spinu 
    +   arg[2]: E_spind
    +   arg[3]: occ_spind
    +   arg[4]: vac (electrostatic potential)
    +   arg[5]: show_values (0: no; 1: yes)
    +   arg[6]: vbm of pristine hBN
    +   arg[7]: cbm of pristine hBN
    ============================================================================
    """
    #plt.axhline(0.0, color="k", linestyle="dashed", linewidth=0.5)
    if len(args) == 1:
        bar_width = 0.5
        pass
    else:
        num_lines = len(args) - 1
        for i in range(num_lines):
            plt.axvline(
                1.0/len(args)*(i+1), color="k", 
                linestyle="solid", linewidth=0.8
                )
    bar_width = 1.0 / len(args) / 2 

    for i, arg in enumerate(args):
        for j, yval in enumerate(arg[0]):
            x = [(i+0.1)*bar_width*2, (i+0.45)*bar_width*2]
            y = [yval-arg[4], yval-arg[4]]
            plt.plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[1][j] == 1:
                    plt.text(
                        (i+0.1)*bar_width*2,
                        yval-arg[4],
                        "\u2b06|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.1)*bar_width*2,
                        yval-arg[4],
                        "\u21e7|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[1][j] == 1:
                    plt.text(
                        (i+0.25+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u2b06",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.25+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u21e7",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        for j, yval in enumerate(arg[2]):
            x = [(i+0.55)*bar_width*2, (i+0.9)*bar_width*2]
            y = [yval-arg[4], yval-arg[4]]
            plt.plot(x, y, linewidth=1, color='black')
            if arg[5] == 1:
                if arg[3][j] == 1:
                    plt.text(
                        (i+0.9)*bar_width*2,
                        yval-arg[4],
                        "\u2b07|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.9)*bar_width*2,
                        yval-arg[4],
                        "\u21e9|" + str((format(yval-arg[4], ".4f"))),
                        fontsize=8,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
            else:
                if arg[3][j] == 1:
                    plt.text(
                        (i+0.7+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u2b07",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
                else:
                    plt.text(
                        (i+0.7+j%2*0.05)*bar_width*2,
                        yval-arg[4],
                        "\u21e9",
                        fontsize=20,
                        horizontalalignment='center', 
                        verticalalignment='center'
                        )
        
        plt.fill_between(
            np.array([i,i+1])*bar_width*2,  # x range to fill
            [arg[7], arg[7]],               # lower limit of y
            #[0, 0],                         # upper limit of y
            [10, 10],                       # upper limit of y
            facecolor="tab:red",             # The fill color
            alpha=0.32,                      # Transparency of the fill
            )
        plt.fill_between(
            np.array([i,i+1])*bar_width*2,  # x range to fill
            [-10, -10],                     # lower limit of y
            [arg[6], arg[6]],               # upper limit of y
            facecolor="tab:blue",               # The fill color
            alpha=0.32                      # Transparency of the fill
            )


def plot_config_coord_diagram(etot_1, etot_2, dQ_1, dQ_2, xlim, ylim, **kwargs):
    """
    ============================================================================
    +   etot_1: the arrays of total energies of state 1
    +   etot_2: the arrays of total energies of state 2
    +   dQ_1: change of nuclear coordinate of state 1
    +   dQ_2: change of nuclear coordinate of state 2
    +   xlim: limitation of x axis for plot
    +   ylim: limitation of y axis for plot
    +   kwargs: arrows, labels
    +   arrows == {"left_arrow_shift": val, "right_arrow_shift": val,
                "elongate": val, "E_zpl_shift": value, "E_rel_shift": val}
    +   labels == {"label1": {"x": xval, "y": yval, "name": "??"},
                "label2": {"x": xval, "y": yval, "name": "??"}}
    ============================================================================
    """
    npoints = len(etot_1)
    min_etot = min(min(etot_1), min(etot_2))
    # set etot_1 and dQ_1 as the energy surface and change of nuclear coordinates of the state with lower energy
    if min(etot_1) - min(etot_2) < 0:
        pass
    else:
        temp_etot = etot_1
        etot_1 = etot_2
        etot_2 = temp_etot
        temp_dQ = dQ_1
        dQ_1 = dQ_2
        dQ_2 = temp_dQ
    sec_min_etot = min(etot_2)
    max_etot = max(etot_2)
    sec_max_etot = max(etot_1)
    for i in range(npoints):
        if etot_1[i] == min_etot:
            x_of_min_etot = dQ_1[i]
        if etot_2[i] == sec_min_etot:
            x_of_sec_min_etot = dQ_2[i]

    E_zpl = float(format(sec_min_etot - min_etot, ".5f"))
    E_rel = float(format(sec_max_etot - min_etot, ".5f"))
    E_abs = float(format(max_etot - min_etot, ".5f"))
    E_em = float(format(sec_min_etot - sec_max_etot, ".5f"))
    sys.stdout.write("\rData from calculation:")
    sys.stdout.write("\rE_zpl = {} eV".format(E_zpl)) # ZPL
    sys.stdout.write("\rE_rel = {} eV".format(E_rel)) # The energy of gs in es geometry
    sys.stdout.write("\rE_abs = {} eV".format(E_abs)) # absorption
    sys.stdout.write("\rE_em = {} eV\n\n".format(E_em)) # emission

    # this part does plotting
    config_plot()

    best_vals_1 = best_vals_of_quadratic_fct(dQ_1[0:5], etot_1[0:5])
    best_vals_2 = best_vals_of_quadratic_fct(dQ_2[7:14], etot_2[7:14])
    min_etot_fix = quadratic_fct(
        x_of_min_etot,
        best_vals_1[0], best_vals_1[1], best_vals_1[2]
        )
    sec_max_etot_fix = quadratic_fct(
        x_of_sec_min_etot,
        best_vals_1[0], best_vals_1[1], best_vals_1[2]
        )
    sec_min_etot_fix = quadratic_fct(
        x_of_sec_min_etot,
        best_vals_2[0], best_vals_2[1], best_vals_2[2]
        )
    max_etot_fix = quadratic_fct(
        x_of_min_etot,
        best_vals_2[0], best_vals_2[1], best_vals_2[2]
        )
    E_rel_fix = float(format(sec_max_etot_fix - min_etot_fix, ".5f"))
    E_zpl_fix = float(format(sec_min_etot_fix - min_etot_fix, ".5f"))

    min_x = min(xlim)
    max_x = max(xlim)
    #min_y = min(ylim)
    #max_y = max(ylim)
    x = np.arange(min_x, max_x, 0.001)
    y_1 = []
    for k in range(len(x)):
        temp_y_1 = quadratic_fct(
            x[k], best_vals_1[0], best_vals_1[1], best_vals_1[2]
            ) - min_etot
        y_1.append(temp_y_1)
    y_2 = []
    for k in range(len(x)):
        temp_y_2 = quadratic_fct(
            x[k], best_vals_2[0], best_vals_2[1], best_vals_2[2]
            ) - min_etot
        y_2.append(temp_y_2)

    if "arrows" in kwargs:
        values = kwargs.get("arrows")
        plt.annotate(
            '', # text
            xy=(
                x_of_min_etot-values["left_arrow_shift"],\
                0-values["elongate"]
                ), # the point (x, y) to annotate
            xytext=(
                x_of_min_etot-values["left_arrow_shift"],\
                E_zpl_fix+values["elongate"]
                ), # the position to place the text at. If None, defaults to (x, y)
            arrowprops=dict(arrowstyle="<|-|>", color = "k")
            )
        plt.hlines(
            E_rel_fix, x_of_sec_min_etot,
            x_of_sec_min_etot+values["right_arrow_shift"],
            linestyles="dashed"
            )
        plt.annotate(
            '',
            xy=(
                x_of_sec_min_etot+values["right_arrow_shift"],\
                0-values["elongate"]
                ),
            xytext=(
                x_of_sec_min_etot+values["right_arrow_shift"],\
                E_rel_fix+values["elongate"]
                ),
            arrowprops=dict(arrowstyle="<|-|>", color = "k")
            )
        plt.text(
            x_of_min_etot-values["E_zpl_shift"],
            E_zpl_fix/2.0, "\u0394E"
            )
        plt.text(
            x_of_sec_min_etot+values["E_rel_shift"],
            E_rel_fix/2.0, "\u0394E$_{rel}$"
            )

    if "labels" in kwargs:
        labels = kwargs.get("labels")
        plt.text(
            labels["label1"]["x"], labels["label1"]["y"],
            labels["label1"]["name"]
            )
        plt.text(
            labels["label2"]["x"], labels["label2"]["y"],
            labels["label2"]["name"]
            )

    plt.plot(x, y_1, linewidth=2, color="tab:blue")
    plt.plot(x, y_2, linewidth=2, color="tab:red")

    for i in range(npoints):
        plt.plot(
            dQ_1[i], etot_1[i]-min_etot, marker="o", markersize=6,
            markerfacecolor="w", color="tab:blue"
            )
        plt.plot(
            dQ_2[i], etot_2[i]-min_etot, marker="o", markersize=6,
            markerfacecolor="w", color="tab:red"
            )

    plt.hlines(0.0, min_x, max_x, linestyles="dashed")
    plt.hlines(E_zpl_fix, min_x, max_x, linestyles="dashed")

    plt.vlines(x_of_min_etot, -10, 10, linestyles="dashed")
    plt.vlines(x_of_sec_min_etot, -10, 10, linestyles="dashed")

    E_abs_fix = float(format(max_etot_fix - min_etot_fix, ".5f"))
    E_em_fix = float(format(sec_min_etot_fix - sec_max_etot_fix, ".5f"))
    sys.stdout.write("\rData from fitting:")
    sys.stdout.write("\rE_zpl_fix = {} eV".format(E_zpl_fix)) # ZPL
    sys.stdout.write("\rE_rel_fix = {} eV".format(E_rel_fix)) # The energy of gs in es geometry
    sys.stdout.write("\rE_abs_fix = {} eV".format(E_abs_fix)) # absorption
    sys.stdout.write("\rE_em_fix = {} eV\n".format(E_em_fix)) # emission

    y_diff = min(np.abs(np.array(y_2)-np.array(y_1)))
    for i in range(len(x)):
        if np.abs(y_2[i] - y_1[i]) == y_diff:
            E_barrier = float(format(y_2[i] - min(y_2), ".6f"))
    # estimated transition rate through energy barrier
    rate_eff = np.power(10.0,12) * np.exp(-E_barrier*ev2J/(kB*T_room))
    time_eff = 1.0/rate_eff
    rate_eff = "{:e}".format(rate_eff, ".5f")
    time_eff = "{:e}".format(time_eff, ".5f")
    sys.stdout.write("\rE_barrier = {} eV".format(E_barrier))
    sys.stdout.write("\rrate_eff = {} s^-1".format(rate_eff))
    sys.stdout.write("\rtime_eff = {} s".format(time_eff))
    sys.stdout.flush()


def view_3d(x, y, z, **kwargs):
    fig = plt.figure(
        num=None, figsize=(10, 7.5), dpi=120,
        facecolor='w', edgecolor='k'
        )
    ax = fig.add_subplot(1,1,1, projection="3d")
    ax.scatter(x, y, z, zdir="z", s=20, c=None)
    surf = ax.plot_trisurf(
        x, y, z, linewidth=0.2,
        antialiased=True, cmap=plt.cm.viridis, alpha=0.6
        )

    if "circle" in kwargs:
        alpha = np.linspace(0.0, np.pi*2, 200)
        (r, x0, y0, z0) = kwargs.get("circle")
        x1 = r * np.cos(alpha) + x0
        y1 = r * np.sin(alpha) + y0
        z1 = 10
        ax.plot(x1, y1, z1)
    else:
        sys.stdout.write(
            "\rTip: to show circle with defects as center, " +
            "specify 'circle = (r, x0, y0, z0)' from input"
            )

    sys.stdout.write("\nPlot \rTips")
    # title
    if "title" in kwargs:
        if kwargs.get("title") == None:
            sys.stdout.write("\rNo title")
        else:
            ax.set_title("{}".format(kwargs.get("title")))
    else:
        sys.stdout.write(
            "\rTip: to add title to diagram, specify " +
            "'title = {}'".format("?") + "from input"
            )

    # about view directions and axes labels
    if "view_direction" in kwargs:
        if kwargs.get("view_direction") == "front_view":
            ax.view_init(azim=-90, elev=0)
            ax.w_yaxis.line.set_lw(0.)
            ax.set_yticks([])
            if "xlabel" in kwargs:
                if kwargs.get("xlabel") == None:
                    sys.stdout.write("\rNo xlabel")
                else:
                    ax.set_xlabel(kwargs.get("xlabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show xlabel, specify 'xlabel = ?' from input"
                    )
            if "zlabel" in kwargs:
                if kwargs.get("zlabel") == None:
                    sys.stdout.write("\rNo zlabel")
                else:
                    ax.set_zlabel(kwargs.get("zlabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show zlabel, specify 'zlabel = ?' from input"
                    )
        elif kwargs.get("view_direction") == "top_view":
            ax.view_init(azim=-90, elev=90)
            ax.w_zaxis.line.set_lw(0.)
            ax.set_zticks([])
            if "xlabel" in kwargs:
                if kwargs.get("xlabel") == None:
                    sys.stdout.write("\rNo xlabel")
                else:
                    ax.set_xlabel(kwargs.get("xlabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show xlabel, specify 'xlabel = ?' from input"
                    )
            if "ylabel" in kwargs:
                if kwargs.get("ylabel") == None:
                    sys.stdout.write("\rNo ylabel")
                else:
                    ax.set_ylabel(kwargs.get("ylabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show ylabel, specify 'ylabel = ?' from input"
                    )
        else:
            ax.view_init(azim=-180, elev=0)
            ax.w_xaxis.line.set_lw(0.)
            ax.set_xticks([])
            if "ylabel" in kwargs:
                if kwargs.get("ylabel") == None:
                    sys.stdout.write("\rNo ylabel")
                else:
                    ax.set_ylabel(kwargs.get("ylabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show ylabel, specify 'ylabel = ?' from input"
                    )
            if "zlabel" in kwargs:
                if kwargs.get("zlabel") == None:
                    sys.stdout.write("\rNo zlabel")
                else:
                    ax.set_zlabel(kwargs.get("zlabel"))
            else:
                sys.stdout.write(
                    "\rTip: to show zlabel, specify 'zlabel = ?' from input"
                    )
    else:
        sys.stdout.write(
            "\rTip: to take the top view of the 3D plot, " +
            "specify 'view_direction = top_view' from input, so as " +
            "for front view and left view."
            )
        if "xlabel" in kwargs:
            if kwargs.get("xlabel") == None:
                sys.stdout.write("\rNo xlabel")
            else:
                ax.set_xlabel(kwargs.get("xlabel"))
        else:
            sys.stdout.write(
                "\rTip: to show xlabel, specify " +
                "'xlabel = ?' from input"
                )
        if "ylabel" in kwargs:
            if kwargs.get("ylabel") == None:
                sys.stdout.write("\rNo ylabel")
            else:
                ax.set_ylabel(kwargs.get("ylabel"))
        else:
            sys.stdout.write(
                "\rTip: to show ylabel, specify 'ylabel = ?' from input"
                )
        if "zlabel" in kwargs:
            if kwargs.get("zlabel") == None:
                sys.stdout.write("\rNo zlabel")
            else:
                ax.set_zlabel(kwargs.get("zlabel"))
        else:
            sys.stdout.write(
                "\rTip: to show zlabel, specify 'zlabel = ?' from input"
                )

    # range of axes
    if "xlim" in kwargs:
        if kwargs.get("xlim") == None:
            sys.stdout.write("\rNo xlim")
        else:
            ax.set_xlim(kwargs.get("xlim"))
    else:
        sys.stdout.write(
            "\rTip: to constrain x-axis, specify 'xlim = (min,max)' from input"
            )
    if "ylim" in kwargs:
        if kwargs.get("ylim") == None:
            sys.stdout.write("\rNo ylim")
        else:
            ax.set_ylim(kwargs.get("ylim"))
    else:
        sys.stdout.write(
            "\rTip: to constrain y-axis, specify 'ylim = (min,max)' from input"
            )
    if "zlim" in kwargs:
        if kwargs.get("zlim") == None:
            sys.stdout.write("\rNo zlim")
        else:
            ax.set_zlim(kwargs.get("zlim"))
    else:
        sys.stdout.write(
            "\rTip: to constrain z-axis, specify 'zlim = (min,max)' from input"
            )

    # colorbar
    if "plot_corlorbar" in kwargs:
        if kwargs.get("plot_corlorbar"):
            if "setup_boundaries" in kwargs:
                fig.colorbar(surf, boundaries=kwargs.get("setup_boundaries"))
            else:
                sys.stdout.write(
                    "\rTip: to set boundaries for colorbar, " +
                    "specify 'setup_boundaries = np.linspace(min,max)' " +
                    "from input"
                    )
                fig.colorbar(surf)
    else:
        sys.stdout.write(
            "\rTip: to plot colorbar, specify " +
            "'plot_corlorbar = True' from input"
            )
    sys.stdout.flush()

