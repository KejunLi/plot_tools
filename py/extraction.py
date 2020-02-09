#!/usr/bin/env python3
import re
import numpy as np
from constants import *


def extract_filename(f_name, condition):
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
    for line in lines:
        for condition in conditions:
            if condition in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line)
                etot = float(raw_etot[0])*Ry2eV
                list_etot.append(etot)
    return(list_etot)



def extract_eps(dir_f):
    """
    this function read eps file with input of directory and file name
    sequentially output the lists of photon energy, Re(eps), Im(eps),
    and directory+file name
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



def extract_cellpara(dir_f):
    """
    Look for celldm from input file, and if there is no celldm in input, look
    for it from output file
    list_cellpara is the CELL_PARAMETERS in cartesian coordinates
    return:
    type(list_cellpara) = list
    list_cellpara = [array1,array2,array3]
    To use list_cellpara, do the following
    [x, y, z] = np.dot(np.asarray(list_cellpara), [crystal coordinates])
    """
    found_CELL_PARAMETERS = False
    got_CELL_PARAMETERS = False
    got_celldm = False
    read_times = int(0)
    with open(dir_f, "r") as f:
        lines = f.readlines()
        list_cellpara = []
        list_celldm = []
    # try read input
    if ".in" in dir_f:
        for line in lines:
            axes = []
            if not found_CELL_PARAMETERS and "CELL_PARAMETERS" not in line:
                continue
            elif not found_CELL_PARAMETERS and "CELL_PARAMETERS" in line:
                if "angstrom" in line:
                    convertunit = 1.0
                elif "bohr" in line:
                    convertunit = Bohr2Ang
                found_CELL_PARAMETERS = True
                continue
            if found_CELL_PARAMETERS and read_times < 3:
                for x in line.strip().split():
                    axes.append(float(x))
                list_cellpara.append(np.dot(np.asarray(axes), convertunit))
                read_times += 1
            else:
                break
    elif ".out" in dir_f:
        for line in lines:
            axes = []
            if not got_celldm and "celldm(1)" not in line:
                continue
            elif not got_celldm and "celldm(1)" in line:
                # unit is Bohr
                std_celldm = float(re.findall(r"[+-]?\d+\.\d*", line)[0])
                got_celldm = True
                continue
            if not found_CELL_PARAMETERS and "crystal" not in line:
                continue
            elif not found_CELL_PARAMETERS and "crystal" in line:
                found_CELL_PARAMETERS = True
                continue
            elif found_CELL_PARAMETERS and read_times < 3:
                for x in re.findall(r"[+-]?\d+\.\d*", line):
                    axes.append(float(x))
                # get CELL_PARAMETERS and convert unit at the same time
                list_cellpara.append(np.dot(np.asarray(axes), \
                    std_celldm*Bohr2Ang))
                read_times += 1
            else:
                break
    return(list_cellpara)



def extract_aps(dir_f):
    """
    this can read relax.in, relax.out or scf.in
    find, read and save the optimized ATOMIC POSITIONS
    output is an array of atomic positions
    l_atom_atompos = [[l_atom], [l_atompos]]
    l_atom = [atom1, atom2, atom3, ...]
    l_atompos = [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], ...]
    """
    found_atompos = False
    got_all_atompos = False
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
        # directly go to the last ATOMIC_POSITIONS
        if counts > 1:
            if "ATOMIC_POSITIONS" in line:
                counts -= 1
                continue
            else:
                continue
        elif counts == 1 and not found_atompos:
            if "ATOMIC_POSITIONS" not in line:
                continue
            elif "ATOMIC_POSITIONS" in line:
                found_atompos = True
                continue
        if found_atompos:
            l_raw_atom_atompos = line.strip().split()
            if len(l_raw_atom_atompos) != 4:
                got_all_atompos = True
            elif not got_all_atompos:
                #print(line)
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


def extract_eigenenergy(dir_f):
    """
    this function is used to extract eigenenergy and occupations at Gamma point
    near valence band maximum (VBM) and conduction band minimum (CBM)
    return:
    l_all = [l_E_spinup, l_E_spindown, l_occ_spinup, l_occ_spindown]
    type(l_E_spinup) = list
    type(l_E_spindown) = list
    type(l_occ_spinup) = list
    type(l_occ_spindown) = list
    """
    found_E_spinup = False
    found_E_spindown = False
    found_kb_spinup = False
    found_kb_spindown = False
    is_occ_spinup = False
    is_occ_spindown = False
    got_E_spinup = int(0)
    got_E_spindown = int(0)
    got_occ_spinup = int(0)
    got_occ_spindown = int(0)
    l_all = []
    l_E_spinup = []
    l_E_spindown = []
    l_occ_spinup = []
    l_occ_spindown = []
    counts = int(0)
    with open(dir_f, "r") as f:
        lines = f.readlines()
    for line in lines:
        if "SPIN UP" in line:
            counts += 1
    num_spinup = counts
    num_spindown = counts
    # print(counts)
    # collect spinup eigenenergies
    for line in lines:
        if num_spinup > 1:
            if "SPIN UP" in line:
                num_spinup -= 1
                continue
            else:
                continue
        elif num_spinup == 1 and not found_E_spinup:
            if "SPIN UP" not in line:
                continue
            elif "SPIN UP" in line:
                found_E_spinup = True
                continue
        elif found_E_spinup:
            if "k = " not in line and not found_kb_spinup:
                continue
            elif "k = " in line and not found_kb_spinup:
                found_kb_spinup = True
                continue
            if found_kb_spinup and not line.split():
                if not got_E_spinup:
                    continue
                elif got_E_spinup and not got_occ_spinup:
                    is_occ_spinup = True
                    continue
                else:
                    break
            elif found_kb_spinup and line.split() and not is_occ_spinup:
                for value in line.strip("\n").split():
                    l_E_spinup.append(float(value))
                got_E_spinup += len(line.strip("\n").split())
            elif is_occ_spinup and line.split():
                if len(line.split()) != 8:
                    continue
                else:
                    for value in line.split():
                        l_occ_spinup.append(float(value))
                    got_occ_spinup += len(line.split())

    # collect spindown eigenenergies
    for line in lines:
        if num_spindown > 1:
            if "SPIN DOWN" in line:
                num_spindown -= 1
                continue
            else:
                continue
        elif num_spindown == 1 and not found_E_spindown:
            if "SPIN DOWN" not in line:
                continue
            elif "SPIN DOWN" in line:
                found_E_spindown = True
                continue
        elif found_E_spindown:
            if "k = " not in line and not found_kb_spindown:
                continue
            elif "k = " in line and not found_kb_spindown:
                found_kb_spindown = True
                continue
            if found_kb_spindown and not line.split():
                if not got_E_spindown:
                    continue
                elif got_E_spindown and not got_occ_spindown:
                    is_occ_spindown = True
                    continue
                else:
                    break
            elif found_kb_spindown and line.split() and not is_occ_spindown:
                for value in line.strip("\n").split():
                    l_E_spindown.append(float(value))
                got_E_spindown += len(line.strip("\n").split())
            elif is_occ_spindown and line.split():
                if len(line.split()) != 8:
                    continue
                else:
                    for value in line.split():
                        l_occ_spindown.append(float(value))
                    got_occ_spindown += len(line.split())

    l_all.append(l_E_spinup)
    l_all.append(l_E_spindown)
    l_all.append(l_occ_spinup)
    l_all.append(l_occ_spindown)
    return(l_all)
