#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from sort_files import files_in_dir
from configuration_for_plot import config_plot
from var_color import var_color
import sys
import os
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
from read_qe import read_vac, qe_out


plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
path = "/home/likejun/work/pristine_hbn/6x6/pbe/vbm_vac"
length = np.arange(7, 61, 1)
vbm = np.zeros(len(length), dtype=float)
vac = np.zeros(len(length), dtype=float)
paths_zlength = []
paths_vac = []
for i in range(len(length)):
    qe = qe_out(
        os.path.join(path, "dzlength_"+str(length[i]), "zlength_"+str(length[i])+".out")
        )
    qe.read_eigenenergies()
    qe.read_bandgap()
    vbm[i] = qe.vbm
    vac[i] = np.amax(
        read_vac(
            os.path.join(path, "dzlength_"+str(length[i]), "job_vacuum/avg.out")
            )[1]
        )
print(vac)
plt.plot(length, vbm, color="tab:red", marker="o", label="$VBM$")
plt.plot(length, vac, color="tab:blue", marker="o", label="$V_{vac}$")
plt.plot(length, vbm-vac, color="black", marker="o", label="$VBM-V_{vac}$")
plt.legend()
plt.xlabel("$L_z$ (length of supercell in z-axis, \u212b)")
plt.ylabel("$E$ (eV)")
plt.title("$VBM$ and $V_{vac}$ vs $L_{z}$")
plt.text(40, -2, "Pristine 6x6 hBN\nnk331, ecutwfc=50 Ry")
plt.show()


