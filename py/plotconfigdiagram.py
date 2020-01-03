#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pylab as ply
import numpy as np
import re
import os
from configuration_for_plot import config_plot
from sort_files import sort_var_and_f

# read *.out to find lines with total energy and extract data from lines
def extract_etot(dir_f):
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


# collect files with satisfactory names in a designated directory
# objective files should contain con_f in their names
def files_in_dir(d_dir, con_f):
    d_dir = d_dir + "/"
    l_f = []
    l_df = []
    for f_name in os.listdir(d_dir):
        if con_f in f_name:
            l_f.append(f_name)
            l_df.append(d_dir+f_name)
    l_f_df = [l_f, l_df]
    #print(l_f)
    #print(l_df)
    return(l_f_df)



directory = "/home/KEJUNLI/work/hBN/Ti/nonradia/my_template_yinan_structure"
dir_lin = []
dir_ratio = []
dir_f = []
dir_lin = files_in_dir(directory, "lin")[1]
for d_lin in dir_lin:
    dir_ratio.append(files_in_dir(d_lin, "ratio-")[1])
#print(dir_ratio)
for d_ratio in dir_ratio:
    dir_f_temp = []
    for d_rat in d_ratio:
        dir_f_temp.append(files_in_dir(d_rat, "scf.out")[1][0])
    dir_f.append(dir_f_temp)
print(dir_f)

for d_f in dir_f:
    l_etot = []
    s_ratio = sort_var_and_f(d_f)[0]
    s_dir_f = sort_var_and_f(d_f)[1]
    print(len(s_ratio), len(s_dir_f))
    for i in range(len(s_dir_f)):
        l_etot.append(extract_etot(s_dir_f[i]))
    color_list = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:pink", "tab:cyan"]
    #plt.plot(s_ratio, l_etot, marker = "", markersize = 2, color = color_list[i])
    plt.plot(s_ratio, l_etot)

#plt.plot([0,1], [l_min_etot[1]-l_min_etot[0],l_min_etot[1]-l_min_etot[0]], "k:")
#plt.text(0.75, 0.2, "E$_{ZPL}$ = "+str(round(l_min_etot[1]-l_min_etot[0], 3))+"eV")
plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title("Kejun")
#plt.xlim([0,1])
plt.show()

"""
config_plot()

a = [0, 0]
l_min_etot = [0, 0]
l_ratio = [0.0000, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.5000, 0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 1.0000]
destination = ["/home/KEJUNLI/work/hBN/Ti/nonradia/my_template_yinan_structure/lin-gs", "/home/KEJUNLI/work/hBN/Ti/nonradia/my_template_yinan_structure/lin-cdftup1"]
plt.plot([0,1], [0,0], "k:")
for i, d in enumerate(destination):
    l_etot = []
    sle = []
    for ratio in l_ratio:
        ratio = format(ratio, ".4f")
        d_dir = d + "ratio-" + str(ratio) + "/" 
        l_etot.append(extract_etot(d_dir, "scf.out"))
    #print(l_ratio)
    #print(l_etot)
    l_min_etot[i] = min(l_etot)[0]
    a[i] = l_etot[12][0]
    print(l_min_etot[i])
    print(l_etot[0][0]-l_min_etot[0])
    print(l_min_etot[1]-a[0])
    for etot in l_etot:
        sle.append(float(etot[0]-l_min_etot[0]))
    color_list = ["tab:red", "tab:orange", "tab:green", "tab:blue", "tab:purple", "tab:pink", "tab:cyan"]
    plt.plot(l_ratio, sle, marker = "", markersize = 2, color = color_list[i+1])
#plt.legend(loc = "upper left")
plt.plot([0,1], [l_min_etot[1]-l_min_etot[0],l_min_etot[1]-l_min_etot[0]], "k:")
plt.text(0.75, 0.2, "E$_{ZPL}$ = "+str(round(l_min_etot[1]-l_min_etot[0], 3))+"eV")
plt.xlabel("\u0394Q (amu$^{1/2}$$\AA$)")
plt.ylabel("E (eV)")
plt.title("Kejun")
plt.xlim([0,1])
plt.show()


"""
