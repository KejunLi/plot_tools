#!/usr/bin/env python3
import re
import numpy as np


def extract_fn(f_name, condition):
    """
    read file name and extract characteristic feature in file name
    output is a list
    list_output = [chara_in_name, num_in_name, chara_name]
    type(chara_in_name) = str
    type(num_in_name) = str
    type(chara_name) = str
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
    output is a list
    list_etot = [etot1, etot2, etot3, ...]
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
    output is an array which contains three lists
    data_set = [E, re_eps, im_eps]
    E = [E1, E2, E3, ...]
    re_eps = [re_eps1, re_eps2, re_eps3, ...]
    im_eps = [im_eps1, im_eps2, im_eps3, ...]
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
    output is a list
    l_cpu_time = [cpu_time1, cpu_time2, cpu_time3, ...]
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


def extract_aps(dir_f):
    """
    this function reads relax.out or scf.out
    find, read and save the optimized ATOMIC POSITIONS
    output is an array of atomic positions
    l_atom_atompos = [[l_atom], [l_atompos]]
    l_atom = [atom1, atom2, atom3, ...]
    l_atompos = [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], ...]
    """
    found_atompos = False
    got_all_atompos = False
    i=0
    with open(dir_f, "r") as f:
        lines = f.readlines()
        l_atom_atompos = []
        counts = int(0)
        l_atom = []
        l_atompos = []
    for line in lines:
        if "ATOMIC_POSITIONS" in line:
            counts += 1
    # print(counts)
    for line in lines:
        l_raw_atompos = []
        atomic_position = []
        if counts > 1:
            if "ATOMIC_POSITIONS" in line:
                counts -= 1
                continue
            else:
                continue
        elif counts == 1 and not found_atompos:
            if "ATOMIC_POSITIONS" not in line:
                i += 1
                # print(found_ap, i)
                continue
            elif "ATOMIC_POSITIONS" in line:
                found_atompos = True
                continue
        if found_atompos:
            l_raw_atom_atompos = line.strip("\n").split()
            if len(l_raw_atom_atompos) != 4:
                got_all_atompos = True
            elif not got_all_atompos:
                # print(line)
                #print(l_raw_atom_atompos)
                atom_name = l_raw_atom_atompos[0]
                x = float(l_raw_atom_atompos[1])
                y = float(l_raw_atom_atompos[2])
                z = float(l_raw_atom_atompos[3])
                atomic_position.append(x)
                atomic_position.append(y)
                atomic_position.append(z)
                l_atompos.append(atomic_position)
                l_atom.append(atom_name)
    l_atom_atompos.append(l_atom)
    l_atom_atompos.append(l_atompos)
    #print(l_atom_atompos)
    return(l_atom_atompos)
