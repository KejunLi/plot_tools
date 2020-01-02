#!/usr/bin/env python3
import os
import re

#####################################################
# d_dir stands for destinated directory
# f_name for name of a file
# l_f for list of files
# l_var for list of variable
# l_dvar for list of dependent variable
# s_var for sorted variable
# sl_var for sorted list of variable
# sl_dvar for sorted list of dependent variable
####################################################

# collect files with satisfactory names in a designated directory
def files_in_dir(d_dir):
    l_f = []
    for f_name in os.listdir(d_dir):
        if "nk" in f_name and ".out" in f_name and "" in f_name:
            l_f.append(f_name)
    #print(l_f)
    return(l_f)


# sort out files with variable in ascending order
# return a sorted list of variable (first column)
# and a sorted list of files (second column)
def sort_var_and_f(l_f):
    l_var = []
    for f_name in l_f:
        raw_var = re.findall(r"[+-]?\d+", f_name)
        var = int(raw_var[0])
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


# refine file name and leave characteristic feature
def refine_f(f_name):
    rf_name = re.findall(r"[q]\d*", f_name)
    # type(rf_name) == <class 'list'>
    return(rf_name[0])







# read *.out to find lines with total energy and extract data from lines
def extract_etot(d_dir, f_name):
    d_file = str(d_dir) + str(f_name)
    with open(d_file, "r") as f:
        lines = f.readlines()
        saved_data = []
        # unit conversion
        ry2ev = 13.6056980659
        for line in lines:
            if "Final" in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line)
                etot = float(raw_etot[0])*ry2ev
                saved_data.append(etot)
    return(saved_data)


# sort out files with variable ascending
def sort_files(d_dir):
    l_var = []
    l_dvar = []
    for f_name in files_in_dir(d_dir):
        raw_var = re.findall(r"[+-]?\d+", f_name)
        var = int(raw_var[0])
        l_var.append(var)
        l_dvar.extend(extract_dvar(d_dir, f_name))
    saved_list = [l_var, l_dvar]

    # create a sorted list of nk and etot
    sl_var = sorted(l_var)
    sl_dvar = []

    # sorted the list of etot according to the order of sl_nk
    for i, s_var in enumerate(sl_var):
        for j, var in enumerate(l_var):
            if var == s_var:
                sl_dvar.append(l_dvar[j])
    sl = [sl_var, sl_dvar]
    print(sl)
    return(sl)

