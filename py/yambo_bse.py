#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from configuration_for_plot import config_plot, plot_eps
from var_color import var_color

config_plot()

plot_eps("/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/pbe0_gwbse_nk221_nbnd1000/data",
"black",
label="IP, nk221, KfnQP_up(dn)_E=2.7(3.6)eV, BSSmod='d'",
#orientation="y",
type="BSE")

plot_eps("/home/likejun/work/tibn/nk331/tibn_oncv_c1/6x6/nonradiative/nbnd1000/data",
"tab:red",
label="BSE, nk331, KfnQP_up(dn)_E=2.7(3.6)eV, BSSmod='h'",
#orientation="y",
type="BSE")



plt.legend(loc = "upper left")
plt.xlabel("E (eV)")
#plt.title(title)
plt.ylabel("Im(${\u03B5}$)")
#plt.xlim(-0.5, 6)
plt.ylim(-0.4, 4.8)
plt.show()
