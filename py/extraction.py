#!/usr/bin/env python3
import re
import numpy as np


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


def extract_aps(dir_f):
    """
    this function reads relax.out or scf.out
    find, read and save the optimized ATOMIC POSITIONS
    output is an array of atomic positions
    l_atom_atompos = [[l_atom], [l_atompos]]
    l_atom = [atom1, atom2, atom3, ...]
    l_atompos = [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], ...]
    """
    isfound_atompos = False
    isgot_all_atompos = False
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
        elif counts == 1 and not isfound_atompos:
            if "ATOMIC_POSITIONS" not in line:
                i += 1
                continue
            elif "ATOMIC_POSITIONS" in line:
                isfound_atompos = True
                continue
        if isfound_atompos:
            l_raw_atom_atompos = line.strip("\n").split()
            if len(l_raw_atom_atompos) != 4:
                isgot_all_atompos = True
            elif not isgot_all_atompos:
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

def extract_eigenenergy():
    """
    this function is used to extract eigenenergy and occupations
    near valence band maximum (VBM) and conduction band minimum (CBM)
    return:
    l_all = [l_E_spinup, l_E_spindown, l_occ_spinup, l_occ_spindown]
    type(l_E_spinup) = array
    type(l_E_spindown) = array
    type(l_occ_spinup) = array
    type(l_occ_spindown) = array
    """
    isfound_E_spinup = False
    isfound_E_spindown = False
    isfound_kb_spinup = False
    isfound_kb_spindown = False
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
        l_raw_E_spinup = []
        if num_spinup > 1:
            if "SPIN UP" in line:
                num_spinup -= 1
                continue
            else:
                continue
        elif num_spinup == 1 and not isfound_E_spinup:
            if "SPIN UP" not in line:
                continue
            elif "SPIN UP" in line:
                isfound_E_spinup = True
                continue
        elif isfound_E_spinup:
            if "k = " not in line and not isfound_kb_spinup:
                continue
            elif "k = " in line and not isfound_kb_spinup:
                isfound_kb_spinup = True
                continue
            if isfound_kb_spinup and not line.split():
                if not got_E_spinup:
                    continue
                elif got_E_spinup and not got_occ_spinup:
                    is_occ_spinup = True
                    continue
                else:
                    break
            elif isfound_kb_spinup and line.split() and not is_occ_spinup:
                l_raw_E_spinup = line.strip("\n").split()
                got_E_spinup += len(l_raw_E_spinup)
                l_E_spinup.append(np.array(l_raw_E_spinup))
            elif is_occ_spinup and line.split():
                if len(line.split()) != 8:
                    continue
                else:
                    l_occ_spinup.append(np.array(line.split()))
                    got_occ_spinup += len(line.split())

    # collect spindown eigenenergies
    for line in lines:
        l_raw_E_spindown = []
        if num_spindown > 1:
            if "SPIN DOWN" in line:
                num_spindown -= 1
                continue
            else:
                continue
        elif num_spindown == 1 and not isfound_E_spindown:
            if "SPIN DOWN" not in line:
                continue
            elif "SPIN DOWN" in line:
                isfound_E_spindown = True
                continue
        elif isfound_E_spindown:
            if "k = " not in line and not isfound_kb_spindown:
                continue
            elif "k = " in line and not isfound_kb_spindown:
                isfound_kb_spindown = True
                continue
            if isfound_kb_spindown and not line.split():
                if not got_E_spindown:
                    continue
                elif got_E_spindown and not got_occ_spindown:
                    is_occ_spindown = True
                    continue
                else:
                    break
            elif isfound_kb_spindown and line.split() and not is_occ_spindown:
                l_raw_E_spindown = line.strip("\n").split()
                got_E_spindown += len(l_raw_E_spindown)
                l_E_spindown.append(np.array(l_raw_E_spindown))
            elif is_occ_spindown and line.split():
                if len(line.split()) != 8:
                    continue
                else:
                    l_occ_spindown.append(np.array(line.split()))
                    got_occ_spindown += len(line.split())
    l_all.append(l_E_spinup)
    l_all.append(l_E_spindown)
    l_all.append(l_occ_spinup)
    l_all.append(l_occ_spindown)
    return(l_all)
