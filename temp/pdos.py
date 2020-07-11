#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from sort_files import files_in_dir
from configuration_for_plot import config_plot
from var_color import var_color

config_plot()


bn = ["b21_s", "b21_p"]
orbit = ["2s", "2p"]
line = ["-.", "-"]
for i in range(len(bn)):
    print(bn[i])
    dir_f = files_in_dir("/home/likejun/work/mobn/mobn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000_qe6.1_yambo4.4/job_pdos", bn[i])[1]
    data = np.loadtxt(dir_f[0], dtype=None)
    E = data[:, 0]
    pdosup = data[:, 1]
    pdosdn = data[:, 2]
    plt.plot(E, pdosdn, color=var_color("tab:blue", 0.3*(i+1)), linestyle=line[i], label=orbit[i] + "(spin down)")
    plt.plot(E, pdosup, color=var_color("tab:red", 0.5*(i+1)), linestyle=line[i], label=orbit[i] + "(spin up)")

plt.legend()
plt.xlim(-9, 9.7)
plt.ylim(0, 2.5)
plt.xlabel("E (eV)")
plt.ylabel("PDOS")
plt.show()

"""
bn = ["b15_s", "b15_p", "ti_4s", "ti_3d"]
orbit = ["B_2s", "B_2p", "Ti_4s", "Ti_3d"]
line = ["--", ":", "-.", "-"]
for i in range(len(bn)):
    print(bn[i])
    dir_f = files_in_dir("/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_exx0.41_gwbse_nk331_nbnd1000_qe6.1_yambo4.4/job_pdos", bn[i])[1]
    data = np.loadtxt(dir_f[0], dtype=None)
    E = data[:, 0]
    pdosup = data[:, 1]
    pdosdn = data[:, 2]
    #plt.plot(E, pdosdn, color=var_color("tab:blue", 0.3*(i+1)), linestyle=line[i], label=orbit[i] + "(spin down)")
    plt.plot(E, pdosup, color=var_color("tab:red", 0.25*(i+1)), linestyle=line[i], label=orbit[i] + "(spin up)")

plt.legend()
plt.xlim(-9, 9.7)
plt.ylim(0, 2.5)
plt.xlabel("E (eV)")
plt.ylabel("PDOS")
plt.show()
"""
