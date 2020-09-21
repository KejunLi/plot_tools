#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, "/home/likejun/work/github/yambo_post_processing")
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from read_gw import read_gw_out
from fit_functions import exponential, best_vals_of_exponential
plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")


path1 = "/home/likejun/work/carbon_dimer/monolayer/6x6/nonradiative"
path2 = "gw_bse_at_pbe_nbnd2000_nk331/convergence/gw_conv_a"
ecut = np.arange(1, 6, 1)
nbnd = np.array([2000])
gap = np.zeros(len(ecut))

for i in range(len(nbnd)):
    for j in range(len(ecut)):
        dir_a = os.path.join(path1, path2, "o-"+str(nbnd[i])+"b_"+str(ecut[j])+"Ry.qp")
        print(dir_a)
        gw = read_gw_out(dir_a)
        gw.read_corr(0)
        print(gw.corr)
        gap[i] = gw.eigenE[0][1] - gw.eigenE[0][0] + gw.corr[0][1] - gw.corr[0][0]


fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)
#ax.plot(x, y, color="tab:red", label="exponential fit")
#ax.text(50, np.mean(gap), r"$\mathrm{y=143339e^{(-0.51x)}-1866}$")
ax.plot(ecut, gap, marker="o", linestyle=':', label="Mo_ONCV_PBE-1.0.upf")
ax.legend()
ax.set_xlabel("ecutwfc (Ry)")
ax.set_ylabel("Energy (eV/atom)")
plt.show()