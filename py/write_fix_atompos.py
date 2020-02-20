#!/usr/bin/env python3
import numpy as np
from configuration_for_plot import view_3d
import matplotlib.pyplot as plt
from general_functions import fix_atompos, read_vasp
from mpl_toolkits.mplot3d import Axes3D

############################## Input ###########################################
dir_f = "/home/likejun/work/hBN/Ti/9x9/nonradiative/relax-gs/relax.in"
dir_f_pristine = "/home/likejun/work/hBN/bulk_hBN_supercell/mp-984_BN/tibn_9x9.vasp"
radius = 0.2
################################################################################
(list_atom, list_atompos, list_atom_coord, list_defect_coord, max_d) = \
fix_atompos(dir_f, radius, "Ti", fix_x=True, fix_y=True)
(list_atompos1, list_atom_coord1) = read_vasp(dir_f_pristine)
number_of_atoms = len(list_atompos)
#print(np.array(list_atompos))
#print(len(list_atompos1))

x = []
y = []
z = []
"""
# before replacement
print("ATOMIC_POSITIONS angstrom")
for i in range(number_of_atoms):
    str_atompos = " "
    for j in range(len(list_atompos[i])):
        str_atompos = str_atompos + format(list_atompos[i][j], ".9f") + "    "
    print(str_atompos)
"""

for i in range(number_of_atoms):
    if len(list_atompos[i]) == 6:
        d = []
        for j in range(number_of_atoms):
            d.append(np.linalg.norm(list_atom_coord1[j]-list_atom_coord[i]))
        min_d = min(d)
        # replace atomic positions of fixed atoms in defect supercell with
        # the atomic positions of the nearest atoms in pristine supercell
        for k in range(number_of_atoms):
            if d[k] == min_d:
                list_atompos[i][0] = list_atompos1[k][0]
                list_atompos[i][1] = list_atompos1[k][1]
                list_atompos[i][2] = list_atompos1[k][2]
    x.append(list_atom_coord[i][0])
    y.append(list_atom_coord[i][1])
    if len(list_atompos[i]) > 3:
        z.append(list_atom_coord[i][2]+10)
    else:
        z.append(list_atom_coord[i][2])

#print(np.array(list_atom_coord))
# after replacement
print("ATOMIC_POSITIONS angstrom")
for i in range(number_of_atoms):
    new_atompos = "  "
    for j in range(len(list_atompos[i])):
        new_atompos += format(list_atompos[i][j], ".9f") + "    "
    atom_atompos = list_atom[i] + new_atompos
    print(atom_atompos)

r = radius * max_d
x0 = list_defect_coord[0][0]
y0 = list_defect_coord[0][1]
z0 = list_defect_coord[0][2]
view_3d(x, y, z, circle=(r, x0, y0, z0), view_direction="top_view")
#view_3d(x, y, z, radius=radius*max_d)

plt.show()
