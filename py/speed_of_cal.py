#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from sort_files import files_in_dir, sort_var_and_f
from extraction import extract_time, extract_filename
from configuration_for_plot import config_plot
from fitting import lin_fct, best_vals_of_lin_fct

#################################### Input ###############################
directory = "/home/likejun/work/dtvsn"
condition = "time"
# type of diagram wanted to plot
#type_of_diagram = "time_vs_iteration"
type_of_diagram = "speed_vs_nodes"
##########################################################################

config_plot()
l_f = files_in_dir(directory, condition)[0]
l_dir_f = files_in_dir(directory, condition)[1]
sl_var = sort_var_and_f(l_f)[0]
sl_f = sort_var_and_f(l_f)[1]
sl_dir_f = []

for i, s_f in enumerate(sl_f):
    for j, f in enumerate(l_f):
        if f == s_f:
            sl_dir_f.append(l_dir_f[j])

# time-iteration diagram
l_speed = []
l_num_nodes = []
for i, dir_f in enumerate(sl_dir_f):
    x = [0,100]
    y = []
    l_cpu_time = sorted(extract_time(dir_f))
    l_iteration = np.arange(1,len(l_cpu_time)+1,1).tolist()
    print(sl_f[i])
    best_vals = best_vals_of_lin_fct(l_iteration, l_cpu_time)
    speed = 1.0/best_vals[1] # speed of calculation
    l_speed.append(speed)
    num_nodes = extract_filename(sl_f[i], "nodes")[1]
    l_num_nodes.append(num_nodes)
    if type_of_diagram == "time_vs_iteration":
        for j, iteration in enumerate(l_iteration):
            plt.plot(iteration, l_cpu_time[j], marker="o", markersize=4)
        for xi in x:
            yi = lin_fct(xi, best_vals[0], best_vals[1])
            y.append(yi)
        plt.plot(x, y, label=extract_filename(sl_f[i], "nodes")[2])
        plt.legend(loc="upper left")
        plt.xlabel("Number of iteration")
        plt.ylabel("CPU time (s)")

#    elif type_of_diagram == "speed_vs_nodes":
#        plt.plot(extract_fn(sl_f[i], "nodes")[1], speed, marker="o")
#        plt.legend(loc="upper left")
#        plt.xlabel("Number of nodes")
#        plt.ylabel("Calculation speed (iterations/s)")
if type_of_diagram == "speed_vs_nodes":
    plt.plot(l_num_nodes, l_speed, color="tab:blue")
    for i, speedi in enumerate(l_speed):
        plt.plot(extract_fn(sl_f[i], "nodes")[1], speedi, marker="o",
                markersize=8, markerfacecolor="white", color="tab:green")
    plt.xlabel("Number of nodes")
    plt.ylabel("Calculation speed (iterations/s)")


#plt.legend(loc="upper left")
#plt.xlabel("Number of iteration")
#plt.ylabel("CPU time (s)")
#plt.xlim([0,70])
#plt.ylim([0,350])
plt.show()
