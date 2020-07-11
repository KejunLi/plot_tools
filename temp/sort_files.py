#!/usr/bin/env python3
import os
import re


def files_in_dir(d_dir, *conditions):
    """
    collect files with satisfactory names in a designated directory
    objective files should satisfy conditions in their names
    """
    d_dir = d_dir + "/"
    list_f = []
    list_dir_f = []
    for f_name in os.listdir(d_dir):
        for condition in conditions:
            if condition in f_name:
                list_f.append(f_name)
                list_dir_f.append(d_dir+f_name)
    sum_list = [list_f, list_dir_f]
    return(sum_list)


def sort_var_and_f(list_f):
    """
    sort out files with variable in ascending order
    return a sorted list of variable (first column)
    and a sorted list of files (second column)
    input is a list of files
    """
    list_var = []
    for f_name in list_f:
        #raw_var = re.findall(r"[+-]?\d+|[+-]?^\d*\.\d*|0", f_name)
        if re.findall(r"\d+\.\d*", f_name):
            raw_var = re.findall(r"\d+\.\d*", f_name)
        elif re.findall(r"\d+", f_name):
            raw_var = re.findall(r"\d+", f_name)
        #print(raw_var)
        var = float(raw_var[-1])
        #print(var)
        list_var.append(var)
    #print(list_var)
    sorted_list_var = sorted(list_var)
    sorted_list_f = []
    for i, s_var in enumerate(sorted_list_var):
        for j, var in enumerate(list_var):
            if var == s_var:
                sorted_list_f.append(list_f[j])
    sum_list = [sorted_list_var, sorted_list_f]
    # print(sl)
    return(sum_list)
