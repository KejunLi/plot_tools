#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from constants import Ry2eV
import sys
sys.path.insert(0, "/home/likejun/work/github/plot_tools")
from plot_tools import plot_ctl
plt.style.use("/home/likejun/work/github/plot_tools/styles/sg-diag")

E_tot = [-13829.815, -13823.986, -13820.721, -13817.110, -13813.549]
E_corr = [7.789, 1.496, 0, 1.313, 5.16]
vac = 0
vbm = -7.8572
cbm = -0.2192
ylim = [-8.2, 0.2]
label = ["-1/-2", "0/-1", "+1/0", "+2/+1"]
title = "Ti-hBN"
ctl_1 = [E_tot, E_corr, vac, vbm, cbm, ylim, label, title]


E_tot = [-14124.219, -14118.485, -14114.953, -14110.633, -14106.453]
E_corr = [7.654, 1.298, 0, 1.287, 5.145]
vac = 0
vbm = -7.8572
cbm = -0.2192
ylim = [-8.4, 1.2]
label = ["-1/-2", "0/-1", "+1/0", "+2/+1"]
title = "Mo-hBN"
ctl = [E_tot, E_corr, vac, vbm, cbm, ylim, label, title]

plot_ctl(ctl_1, ctl)

plt.show()