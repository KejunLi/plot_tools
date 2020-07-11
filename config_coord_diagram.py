#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import sys
from configuration_for_plot import plot_config_coord_diagram
from sort_files import sort_var_and_f
from extraction import extract_etot
from general_functions import get_dQ_from_scf

################################### Input ######################################
directory = "/home/likejun/work/nbvn/nonradiative"
xlim = (-1.2, 1.8)
ylim = (-0.1, 3.8)
arrows = {"left_arrow_shift": 0.15, "right_arrow_shift":0.2,
        "elongate": 0.005, "E_zpl_shift": 0.3, "E_rel_shift": 0.25}
labels = {"label1": {"x": 0.8, "y": 0.05, "name": "singlet ES Cs 2"},
        "label2": {"x": -0.8, "y": 0.3, "name": "triplet ES Cs"}
        }
title = ""
# 1: no labels, no arrow; 2: labels, no arrows
# 3: no labels, arrows; 4: labels, arrows
style = 3
###############################################################################

(set_dir_scfout, dQ) = get_dQ_from_scf(directory)
# this part refines the scf.out files and extracts
# the ratio of linear extrapolation and corresponding total energies
set_etot = []
set_dQ = [] # nuclear coordinate
for i, list_dir_scfout in enumerate(set_dir_scfout):
    #sys.stdout.write(d_f)
    list_etot = []
    nuc_coord = [] # nuclear coordinate
    sorted_list_ratio = sort_var_and_f(list_dir_scfout)[0]
    sorted_list_dir_scfout = sort_var_and_f(list_dir_scfout)[1]
    for j, dir_f in enumerate(sorted_list_dir_scfout):
        nuc_coord.append(sorted_list_ratio[j]*dQ)
        list_etot.append(extract_etot(dir_f, "!")[0])
    set_dQ.append(nuc_coord)
    set_etot.append(list_etot)

if style == 1:
    plot_config_coord_diagram(set_etot[0], set_etot[1], set_dQ[0], set_dQ[1],
        xlim, ylim)
elif style == 2:
    plot_config_coord_diagram(set_etot[0], set_etot[1], set_dQ[0], set_dQ[1],
        xlim, ylim, labels=labels)
elif style == 3:
    plot_config_coord_diagram(set_etot[0], set_etot[1], set_dQ[0], set_dQ[1],
        xlim, ylim, arrows=arrows)
else:
    plot_config_coord_diagram(set_etot[0], set_etot[1], set_dQ[0], set_dQ[1],
        xlim, ylim, arrows=arrows, labels=labels)


plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title(title)
plt.xlim(xlim)
plt.ylim(ylim)
plt.show()
