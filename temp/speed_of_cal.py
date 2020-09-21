#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, "/home/likejun/work/github/qe_post_processing")
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from read_qe import qe_out
plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")

home = "/home/likejun/work/nbvn_vb/sg15_oncv/out-of-plane/nbvn"
cell = ["6x6", "7x7", "8x8", "9x9"]
cells = np.arange(6, 10, 1)
zpl = np.zeros(len(cells))

for i in range(len(cells)):
    dir_gs = os.path.join(home, cell[i], "nonradiative", "relax-gs/relax.out")
    dir_es = os.path.join(home, cell[i], "nonradiative", "relax-cdftup1/relax.out")
    qe_gs = qe_out(dir_gs)
    qe_es = qe_out(dir_es, show_details=True)
    qe_gs.read_etot()
    qe_es.read_etot()
    zpl[i] = - np.amin(qe_gs.etot) + np.amin(qe_es.etot)
print(zpl)
plt.plot(cells, zpl, marker="o")
plt.show()
