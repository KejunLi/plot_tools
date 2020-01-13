#!/usr/bin/env python3
import re
import numpy as np


def extract_fn(f_name, condition):
    """
    read file name and extract characteristic feature in file name
    """
    num_in_name = re.findall(r"\d+", f_name)[0]
    chara_in_name = re.findall(r"{}".format(condition), f_name)[0]
    # type(*_in_name) == <class 'list'>
    chara_name = chara_in_name + " " + num_in_name
    list_output = [chara_in_name, num_in_name, chara_name]
    return(list_output)



def extract_etot(dir_f, *conditions):
    """
    read *.out to find lines with total energy and extract data from lines
    conditions can be "!", "!!" and "Final"
    """
    with open(dir_f, "r") as f:
        lines = f.readlines()
    list_etot = []
    # unit conversion
    ry2ev = 13.6056980659
    for line in lines:
        for condition in conditions:
            if condition in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line) 
                etot = float(raw_etot[0])*ry2ev
                list_etot.append(etot)
    return(list_etot)



def extract_eps(dir_f):
    """
    this function read eps file with input of directory and file name
    sequentially output the lists of photon energy, Re(eps), Im(eps), and directory+file name
    """
    raw_data = []
    E = []
    re_eps = []
    im_eps = []
    data_set = []
    E = np.loadtxt(dir_f, usecols=(0,))
    re_eps = np.loadtxt(dir_f, usecols=(2,))
    im_eps = np.loadtxt(dir_f, usecols=(1,))
    data_set = [E, re_eps, im_eps]
    return(data_set)

def extract_time(dir_f):
    """
    this function reads the total cpu time spent and saves cpu time as a list
    plot the list versus the number of iteration
    the gradient of the curve suggests how fast a calculation is
    """
    with open(dir_f, "r") as f:
        lines = f.readlines()
        l_cpu_time = []
    #for line in lines[::-1]:
    for line in reversed(lines):
        if "total cpu time" in line:
            cpu_time = re.findall(r"[+-]?\d+\.\d*", line)
            l_cpu_time.append(float(cpu_time[0]))
    return(l_cpu_time)
