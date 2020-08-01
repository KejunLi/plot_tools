#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from constants import Ry2eV
from configuration_for_plot import config_plot

dir_ti_y = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd600_qe6.1_yambo4.4/data-20ry/o-bse.alpha_q1_haydock_bse"

config_plot()

data_bn = np.genfromtxt(dir_ti_y, dtype=float)
E_bn = data_bn[:, 0]
eps_bn = data_bn[:, 1]
eps0_bn = data_bn[:, 3]

plt.legend()
plt.plot(E_bn, eps_bn)
plt.show()