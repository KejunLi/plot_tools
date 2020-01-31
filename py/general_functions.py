#1/usr/bin/env/ python3
import numpy as np
import re


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
def cal_dQ(dx, dy, dz, element):
    """
    this function calculates the change of nuclear coordinate
    input is the change of nuclear cartesian coordinate in x-, y-, z-axes and
    species of element name
    type(dQ) = float
    """
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
