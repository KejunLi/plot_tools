#!/usr/bin/env python3
import numpy as np
from configuration_for_plot import view_3d
import matplotlib.pyplot as plt
from general_functions import fix_atompos

dir_f = "/home/likejun/work/hBN/Ti/9x9/nonradiative/relax-gs/relax.out"

(lap, list_atom_coord) = fix_atompos(dir_f, 0.45, "Ti",  fix_y=True)

x = []
y = []
z = []
for i in range(len(list_atom_coord)):
    x.append(list_atom_coord[i][0])
    y.append(list_atom_coord[i][1])
    if len(lap[i]) > 3:
        z.append(list_atom_coord[i][2]+10)
    else:
        z.append(list_atom_coord[i][2])
view_3d(x,y,z)
plt.show()
