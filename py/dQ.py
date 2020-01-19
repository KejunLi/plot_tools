#!/usr/bin/env python3
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from configuration_for_plot import config_plot
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_aps

################################### Input ######################################
directory = "/home/likejun/work/hBN/Ti/supercell_66/nonradia/my_template_my_structure"
filename = "scf.in"
################################################################################

# this part looks for all the scf.out files and save in the list
# for ground state and excited state, respectively.
l_dir_lin = []
l_dir_ratio = []
l_dir_f = []
set_dir_f = []
l_dir_lin = files_in_dir(directory, "lin")[1]
for dir_lin in l_dir_lin:
    l_dir_ratio.append(files_in_dir(dir_lin, "ratio-")[1])
# print(dir_ratio)
for dir_ratio in l_dir_ratio:
    l_dir_f_temp = []
    for dir_ratio_i in dir_ratio:
        l_dir_f_temp.append(files_in_dir(dir_ratio_i, filename)[1][0])
    set_dir_f.append(l_dir_f_temp)

ith_set_atompos = []
ith_set_atom = []
set_atompos = []
set_atom = []
for i, l_dir_f in enumerate(set_dir_f):
    sl_ratio = sort_var_and_f(l_dir_f)[0]
    sl_dir_f = sort_var_and_f(l_dir_f)[1]
    # print(s_dir_f)
    for j in range(len(sl_dir_f)):
        #l_atom = []
        #l_atompos = []
        #print(s_dir_f[j])
        l_atom = extract_aps(sl_dir_f[j])[0]
        l_atompos = extract_aps(sl_dir_f[j])[1]
        ith_set_atom.append(l_atom)
        ith_set_atompos.append(l_atompos)
    set_atom.append(ith_set_atom)
    set_atompos.append(ith_set_atompos)
print(set_atompos)

for i in range(len(set_atompos)):
    if "gs" in set_dir_f[i][0]:
        for j in range(len(set_atompos[i])):
            
