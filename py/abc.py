#!/usr/bin/env python3
import matplotlib.pyplot as plt
from configuration_for_plot import config_plot

x = [1, 2, 3, 4, 5, 6, 7, 8]
y = [0.497442,
0.497541,
0.497566,
0.497582,
0.497552,
0.497543,
0.497528,
0.497542]
config_plot()

plt.plot(x,y,marker='o',markersize=5,linewidth=2)
plt.xlabel('nk nk 1')
plt.ylabel('ZPL (eV)')
plt.show()

