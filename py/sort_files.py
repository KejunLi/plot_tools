#!/usr/bin/env python3
import os
import re

#####################################################
# d_dir stands for destinated directory
# f_name for name of a file
# l_f for list of files
# l_df for list of directory + file name
# l_var for list of variable
# l_dvar for list of dependent variable
# s_var for sorted variable
# sl_var for sorted list of variable
# sl_dvar for sorted list of dependent variable
####################################################

def files_in_dir(d_dir, *con_fs):
    """
    collect files with satisfactory names in a designated directory
    objective files should contain con_f in their names
    """
    d_dir = d_dir + "/"
    l_f = []
    l_df = []
    for f_name in os.listdir(d_dir):
        for con_f in con_fs:
            if con_f in f_name:
                l_f.append(f_name)
                l_df.append(d_dir+f_name)
    l_f_df = [l_f, l_df]
    #print(l_f)
    #print(l_df)
    return(l_f_df)


def sort_var_and_f(l_f):
    """
    sort out files with variable in ascending order
    return a sorted list of variable (first column)
    and a sorted list of files (second column)
    input is a list of files
    """
    l_var = []
    for f_name in l_f:
        #raw_var = re.findall(r"[+-]?\d+|[+-]?^\d*\.\d*|0", f_name)
        if re.findall(r"\d+\.\d*", f_name):
            raw_var = re.findall(r"\d+\.\d*", f_name)
        elif re.findall(r"\d+", f_name):
            raw_var = re.findall(r"\d+", f_name)
        #print(raw_var)
        var = float(raw_var[0])
        #print(var)
        l_var.append(var)
    #print(l_var)
    sl_var = sorted(l_var)
    sl_f = []
    for i, s_var in enumerate(sl_var):
        for j, var in enumerate(l_var):
            if var == s_var:
                sl_f.append(l_f[j])
    sl = [sl_var, sl_f]
    # print(sl)
    return(sl)



