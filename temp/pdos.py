#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from fit_functions import lin_fct, best_vals_of_lin_fct
import sys
import os
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
from read_qe import qe_out


plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
path = "/home/likejun/work/pristine_hbn/6x6/pbe/vbm_vac"
Lz = np.arange(7, 61, 1) # 7 angstroms <= Lz <= 60 angstroms
time = np.zeros(len(Lz), dtype=float)
dense_grid = np.zeros(len(Lz), dtype=float)
fft = np.zeros((len(Lz), 3), dtype=float)

for i in range(len(Lz)):
    qe = qe_out(
        os.path.join(path, "dzlength_"+str(Lz[i]), "zlength_"+str(Lz[i])+".out")
        )
    qe.read_miscellus()
    time[i] = qe.wall_time
    fft[i, :] = qe.fft
    dense_grid[i] = qe.dense_grid

c0, c1 = best_vals_of_lin_fct(Lz, time)
b0, b1 = best_vals_of_lin_fct(Lz, fft[:, 2])
x = np.arange(7, 60)
ytime = c0 + c1*x
yfft = b0 + b1*x
plt.plot(x, ytime, color="k", linestyle="--")
plt.plot(x, yfft, color="k", linestyle="--")

plt.plot(Lz, time, color="tab:red", marker="o", label="wall_time (s)")
plt.plot(Lz, fft[:,2], color="tab:blue", marker="o", label="FFT_z")


plt.legend()
plt.xlabel("$L_z$ (length of supercell in z-axis, \u212b)")
plt.ylabel("")
plt.text(40, 99, "$t_{wall} = (-10.4395 + 3.9917 * L_z) s$")
plt.text(16, 350, "$FFT_{z} = 2.4182 + 8.6592 * L_z$")
plt.show()


