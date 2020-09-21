#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from read_qe import qe_out
plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")



home = "/home/likejun/work/nbvn_vb/sg15_oncv/out-of-plane/nbvn/ecut"
#cell = ["6x6", "7x7", "8x8", "9x9"]
ecut = np.arange(30, 95, 5)
#cells = np.arange(6, 10, 1)
zpl = np.zeros(len(ecut))



for i in range(len(ecut)):
    dir_gs = os.path.join(home, "gs", "decut_"+str(ecut[i]), "ecut_"+str(ecut[i])+".out")
    dir_es = os.path.join(home, "es", "decut_"+str(ecut[i]), "ecut_"+str(ecut[i])+".out")
    qe_gs = qe_out(dir_gs)
    qe_es = qe_out(dir_es)
    qe_gs.read_etot()
    qe_es.read_etot()
    zpl[i] = - np.amin(qe_gs.etot) + np.amin(qe_es.etot)


plt.axhline(2.4107, linestyle="dashed", linewidth=1, color="k")
plt.plot(ecut, zpl, marker="o", label="6x6x1 $\mathrm{N_BV_N}$\n3x3x1 k-point mesh")

#plt.axhline(0.42368, linestyle="dashed", linewidth=1, color="k")
#plt.plot(ecut, zpl, marker="o", label="6x6x1 $\mathrm{V_B}$")
print(zpl)
plt.ylabel("ZPL (eV)")
plt.xlabel("ecutwfc (Ry)")
plt.legend()
plt.title("3x3x1 k-point mesh, SG15 ONCV")
plt.show()