#!/usr/bin/env python3
import re


def extract_fn(f_name):
    """
    read file name and extract characteristic feature in file name
    """
    rf_name = re.findall(r"[q]\d*", f_name)
    # type(rf_name) == <class 'list'>
    return(rf_name[0])



def extract_etot(dir_f):
    """
    read *.out to find lines with total energy and extract data from lines
    """
    with open(dir_f, "r") as f:
        lines = f.readlines()
        saved_data = []
        # unit conversion
        ry2ev = 13.6056980659
        for line in lines:
            if "!" in line:
                # \d +  # the integral part
                # \.    # the decimal point
                # \d *  # some fractional digits
                raw_etot = re.findall(r"[+-]?\d+\.\d*", line) 
                etot = float(raw_etot[0])*ry2ev
                saved_data.append(etot)
    return(saved_data)


