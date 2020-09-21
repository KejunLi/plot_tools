#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
from read_qe import qe_out
import os

plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
dir1 = "/home/likejun/work/nbvn/out-of-plane/ecut"

ecut = np.arange(30, 95, 5)

etot = np.zeros(len(ecut))

for i in range(len(ecut)):
    path_gs = os.path.join(dir1, "gs", "decut_"+ str(ecut[i]), "ecut_"+str(ecut[i])+".out")
    path_es = os.path.join(dir1, "es", "decut_"+ str(ecut[i]), "ecut_"+str(ecut[i])+".out")
    qe_gs = qe_out(path_gs, show_details=True)
    qe_es = qe_out(path_es, show_details=True)
    qe_gs.read_etot()
    qe_es.read_etot()
    etot[i] = qe_es.etot[-1] - qe_gs.etot[-1]

plt.axhline(0.4236, linestyle="dashed", linewidth=1, color="k")
plt.plot(ecut, etot, marker="o", label="$\mathrm{V_B}$")

plt.legend()
#plt.legend(loc = "upper left")
plt.xlabel("ecutwfc (Ry)")
#plt.title(title)
plt.ylabel("ZPL (eV)")

plt.show()
