#!/usr/bin/env python

from plotbands import bndplot
import matplotlib.pyplot as plt


datafile = "bands.dat.gnu"  # bandx.dat.gnu 
fermi = 6.2564              # eV
symmetryfile = "bands.out"  # bandx.out
subplot = plt.subplot(111)  # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.subplot.html
label = ""

bndplot(datafile,fermi,symmetryfile,subplot,label)
plt.show()
