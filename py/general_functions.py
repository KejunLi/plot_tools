#1/usr/bin/env/ python3
import numpy as np
import re
from extraction import extract_aps, extract_cellpara
from sort_files import files_in_dir
import sys


def interpolate(inp, new_length):
    length = len(inp)
    ratio = (length - 1.0) / (new_length - 1.0)
    new_l_data = []
    for i in range(new_length):
        int_part = int(i*ratio) // 1
        fractional_part = int(i*ratio) % 1
        if fractional_part > 0:
            j = int_part + 1
        else:
            j = int_part
        new_data = (1.0 - fractional_part) * inp[int_part] + \
                fractional_part * inp[j]
        new_l_data.append(new_data)
    return(np.array(new_l_data))



def cal_dQ(dx, dy, dz, element):
    """
    this function calculates the change of nuclear coordinate
    input is the change of nuclear cartesian coordinate in x-, y-, z-axes and
    species of element name
    type(dQ) = float
    """
    dict_atom_mass = {"H": 1.008, "He": 4.003, "Li": 6.94, "Be": 6.9012,
                    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
                    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
                    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
                    "Cl": 35.45, "Ar": 39.95, "K": 39.098, "Ca": 40.078,
                    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
                    "Mn": 54.938, "Fe": 55.845, "Co": 58.993, "Ni": 58.693,
                    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630,
                    "As": 74.9922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
                    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
                    "Nb": 101.07, "Mo": 95.95, "Tc": 97, "Ru": 101.91,
                    "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
                    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60,
                    "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33,
                    "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
                    "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
                    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26,
                    "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49,
                    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23,
                    "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
                    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209,
                    "At": 210, "Rn": 222}
    if element in dict_atom_mass:
        dQ = np.sqrt(np.power(dx,2)+np.power(dy,2)+np.power(dz,2)) * \
                    np.sqrt(dict_atom_mass.get(element))
    return(dQ)


def crystal_coord_to_cartesian_coord(a, b, c, u, v, w):
    """
    convert fractional crystal coordinates to cartesian coordinates
    type(xyz) = array
    """
    x = u * a + v * b * np.cos(np.pi*2.0/3.0)
    y = v * b * np.sin(np.pi*2.0/3.0)
    z = w * c
    xyz = np.array([x, y, z])
    return(xyz)


def energy_level(list_E_spinu, list_E_spind, list_occ_spinu, list_occ_spind,
        lowerlimit, upperlimit):
    """
    post process the extracted energy levels of both spin up and spin down
    find the vbm, cbm and corresponding level numbers
    obtain energy levels within the range of [-lowerlimit, upperlimit] in the
    vinicity of vbm

    the input is the list of spinup energies, spindown energies, occupations of
    spinup, occupations of spindown and range of energy levels
    """
    for i, occupation in enumerate(list_occ_spinu):
        if occupation == 0:
            HO_spinu = i - 1
            LU_spinu = i
            break
    for i, occupation in enumerate(list_occ_spind):
        if occupation == 0:
            HO_spind = i - 1
            LU_spind = i
            break

    sys.stdout.write("Plot energy levels from {} to {}."\
        .format(HO_spinu+1-lowerlimit, HO_spinu+upperlimit))
    # find vbm, cbm and fermi level
    if list_E_spinu[HO_spinu] >= list_E_spind[HO_spind]:
        vbm = list_E_spinu[HO_spinu]
        sys.stdout.write("\rspin up {}: the highest occupied state {} eV"\
        .format(HO_spinu+1, vbm))
    else:
        vbm = list_E_spind[HO_spind]
        sys.stdout.write("\rspin down {}: the highest occupied state {} eV"\
        .format(HO_spind+1, vbm))

    if list_E_spinu[LU_spinu] <= list_E_spind[LU_spind]:
        cbm = list_E_spinu[LU_spinu]
        sys.stdout.write("\rspin up {}: the lowest unoccupied state {} eV"\
        .format(LU_spinu+1, cbm))
    else:
        cbm = list_E_spind[LU_spind]
        sys.stdout.write("\rspin down {}: the lowest unoccupied state {} eV"\
        .format(LU_spind+1, cbm))
    sys.stdout.flush()
    # obtain the range of energy levels
    E_spinu = list_E_spinu[range(HO_spinu-lowerlimit, HO_spinu+upperlimit)]
    E_spind = list_E_spind[range(HO_spinu-lowerlimit, HO_spinu+upperlimit)]
    occ_spinu = list_occ_spinu[range(HO_spinu-lowerlimit, HO_spinu+upperlimit)]
    occ_spind = list_occ_spind[range(HO_spinu-lowerlimit, HO_spinu+upperlimit)]
    bands_range = np.array(range(HO_spinu+1-lowerlimit, HO_spinu+1+upperlimit))
    return(E_spinu, E_spind, occ_spinu, occ_spind, bands_range, vbm, cbm)


def fix_atompos(dir_f, radius, defect, **kwargs):
    """
    this function works for 1D, 2D and 3D
    set point defect (e.g. Ti) as center of a circle with radius as input
    measure distance between pristine atoms (e.g. B or N) and defect
    for atoms out of the circle, append "0 0 1" (x_fixed, y_fixed, z_free)
    after fractional crystal coordinates so that the positions of those atoms
    will be fixed, and those in the circle will still be free
    return:
    (list_atom, new_list_atompos, list_atom_coord, list_defect_coord, max_d)
    type(list_atom) = list
    new_list_atompos = [[atompos1], [atompos2], [atompos3], ...]
    type(atompos1) = list
    list_atom_coord = [[atomcoord1], [atomcoord2], [atomcoord3], ...]
    type(atomcoord1) = list
    list_defect_coord = [[defect_coord1], [defect_coord2], [defect_coord3], ...]
    type(max_d) = float
    """
    sys.stdout.write("Tip: to fix atomic positions of atoms " +
        "that are out of the circle, which is centered at defects and has " +
        "radius {}, ".format(radius) + "specify 'fix_pos = fix_x " +
        "(or fix_y, fix_z)' from the input")
    list_cellpara = extract_cellpara(dir_f)
    list_atom = extract_aps(dir_f)[0]
    list_atompos = extract_aps(dir_f)[1]
    new_list_atompos = []
    temp_new_list_atompos = []
    list_defect_num = []
    list_defect_coord = []
    list_atom_coord = []
    list_distance = []
    for i, atom in enumerate(list_atom):
        if defect == atom:
            list_defect_num.append(i)
            defect_coord = np.matmul(list_atompos[i], list_cellpara)
            list_defect_coord.append(defect_coord)
    for atompos in list_atompos:
        (x, y, z) = np.matmul(atompos, list_cellpara)
        list_atom_coord.append(np.asarray([x, y, z]))
    number_of_defects = len(list_defect_num)
    number_of_atoms = len(list_atompos)

    # look for distance between pristine atoms and defects
    for i in list_defect_num:
        temp_list_distance = []
        for j in range(len(list_atom_coord)):
            d = np.linalg.norm(list_atom_coord[j]-list_atom_coord[i])
            z2 = np.power(list_atom_coord[j][2]-list_atom_coord[i][2], 2.0)
            d = np.sqrt(np.power(d, 2.0)-z2)
            temp_list_distance.append(d)
        list_distance.append(temp_list_distance)
    max_d = max(np.array(list_distance).flatten())

    # fix atoms out of circle
    for i in range(number_of_defects):
        temp_list_atompos = []
        for j in range(number_of_atoms):
            d = list_distance[i][j]
            norm_d = d / max_d
            if norm_d > radius:
                if "fix_x" in kwargs:
                    if kwargs.get("fix_x") == True:
                        list_atompos[j].append(0)
                    else:
                        list_atompos[j].append(1)
                else:
                    list_atompos[j].append(1)
                if "fix_y" in kwargs:
                    if kwargs.get("fix_y") == True:
                        list_atompos[j].append(0)
                    else:
                        list_atompos[j].append(1)
                else:
                    list_atompos[j].append(1)
                if "fix_z" in kwargs:
                    if kwargs.get("fix_z") == True:
                        list_atompos[j].append(0)
                    else:
                        list_atompos[j].append(1)
                else:
                    list_atompos[j].append(1)
            temp_list_atompos.append(list_atompos[j])
        temp_new_list_atompos.append(temp_list_atompos)

    if number_of_defects > 1:
        for j in range(number_of_atoms):
            x = max(len(temp_new_list_atompos[0:][j]))
            for i in range(number_of_defects):
                if len(temp_new_list_atompos[i][j]) == x:
                    new_list_atompos.append(temp_new_list_atompos[i][j])
                    continue
    else:
        new_list_atompos = temp_new_list_atompos[0]
    return(list_atom, new_list_atompos, list_atom_coord, list_defect_coord, \
        max_d)



def read_vasp(dir_f):
    """
    read CELL_PARAMETERS and ATOMIC_POSITIONS in vasp file
    return:
    (list_atompos, list_atom_coord)
    """
    list_cellpara = []
    list_atompos = []
    list_atom_coord = []
    with open(dir_f, "r") as f:
        lines = f.readlines()
    for line in lines[2:5]:
        list_cellpara_component = []
        for x in re.findall(r"[+-]?\d+\.\d*", line):
            list_cellpara_component.append(float(x))
        list_cellpara.append(list_cellpara_component)
    for line in lines[8:]:
        atomic_position = []
        u = float(line.strip().split()[0])
        v = float(line.strip().split()[1])
        w = float(line.strip().split()[2])
        atomic_position.append(u)
        atomic_position.append(v)
        atomic_position.append(w)
        list_atompos.append(atomic_position)
    for i in range(len(list_atompos)):
        (x, y, z) = np.matmul(list_atompos[i], list_cellpara)
        list_atom_coord.append(np.asarray([x, y, z]))
    return(list_atompos, list_atom_coord)


def get_dQ_from_scf(directory):
    # this part looks for all the scf.out files and save in the list
    # for ground state and excited state, respectively.
    list_dir_lin = []
    list_dir_ratio = []
    set_dir_scfin = []
    set_dir_scfout = []
    list_dir_lin = files_in_dir(directory, "lin")[1]
    for dir_lin in list_dir_lin:
        list_dir_ratio.append(files_in_dir(dir_lin, "ratio-")[1])
    # sys.stdout.write(dir_ratio)
    for dir_ratio in list_dir_ratio:
        list_dir_scfout_temp = []
        list_dir_scfin_temp = []
        for dir_ratio_i in dir_ratio:
            list_dir_scfin_temp.append(
                        files_in_dir(dir_ratio_i, "scf.in")[1][0])
            list_dir_scfout_temp.append(
                        files_in_dir(dir_ratio_i, "scf.out")[1][0])
        set_dir_scfin.append(list_dir_scfin_temp)
        set_dir_scfout.append(list_dir_scfout_temp)
    # sys.stdout.write(dir_f)

    # calculate dQ
    set_atom = []
    set_atompos = []
    # look for CELL_PARAMETERS from output
    for dir_f in set_dir_scfout[0]:
        if "ratio-0.0000" in dir_f or "ratio-1.0000" in dir_f:
            list_cellpara = np.asarray(extract_cellpara(dir_f))
    # look for atomic positions from scf.in
    for dir_f in set_dir_scfin[0]:
        if "ratio-0.0000" in dir_f or "ratio-1.0000" in dir_f:
            list_atom = extract_aps(dir_f)[0]
            list_atompos = extract_aps(dir_f)[1]
            set_atom.append(list_atom)
            set_atompos.append(np.asarray(list_atompos))
    set_dQ2 = []
    for i in range(len(set_atompos)):
        for j in range(len(set_atompos[i])):
            (x0, y0, z0) = np.matmul(set_atompos[i][j], list_cellpara)
            (xi, yi, zi) = np.matmul(set_atompos[i+1][j], list_cellpara)
            dQi = cal_dQ(xi-x0, yi-y0, zi-z0, set_atom[i][j])
            set_dQ2.append(dQi**2)
        break
    dQ = np.sqrt(sum(set_dQ2))
    sys.stdout.write("ΔQ = {}".format(dQ))
    sys.stdout.flush()
    return(set_dir_scfout, dQ)
