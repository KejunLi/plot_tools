#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from constants import Ry2eV
from configuration_for_plot import config_plot

def neg(e, x, vbm, fac):
    y = e - (np.array(x) + vbm) * fac
    return(y)

def pos(e, x, vbm, fac):
    y = e + (np.array(x) + vbm) * fac
    return(y)


config_plot()
# 6x6 TiBN
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781 - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 #pbe0
#cbm = 1.9212 # pbe0
#e0 = -1015.8193557141 * Ry2eV
#e_n1 = -1015.8388465050 * Ry2eV + (-0.06223571) * Ry2eV * 2
#e_p1 = -1015.6860995243 * Ry2eV + (0.11072122) * Ry2eV * 2
#e_n2 = -1015.7876453456 * Ry2eV + (-0.07669555) * Ry2eV * 2
#e_p2 = -1015.4371388288 * Ry2eV + (0.26940724) * Ry2eV * 2

# 7x7 TiBN
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781 - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 #pbe0
#cbm = 1.9212 # pbe0
#e0 = -1350.4154160489 * Ry2eV
#e_n1 = -1350.4296751426 * Ry2eV + (-0.06267151) * Ry2eV * 2
#e_p1 = -1350.2823053798 * Ry2eV + (0.11136076) * Ry2eV * 2
#e_n2 = -1350.3821428923 * Ry2eV + (-0.07810308) * Ry2eV * 2
#e_p2 = -1350.0364486179 * Ry2eV + (0.27119277) * Ry2eV * 2

# 6x6 MoBN
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781 - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 #pbe0
#cbm = 1.9212 # pbe0
#e0 = -1037.4433783614 * Ry2eV
#e_n1 = -1037.4878936465 * Ry2eV + (-0.06403590) * Ry2eV * 2
#e_p1 = -1037.2621736520 * Ry2eV + (0.11145359) * Ry2eV * 2
#e_n2 = -1037.4364202609 * Ry2eV + (-0.07846342) * Ry2eV * 2
#e_p2 = -1036.9890571664 * Ry2eV + (0.27122071) * Ry2eV * 2

# 6x6 CBVN in-plane qe+jdftx
vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
cbm = 1.1781  - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 # pbe0
#cbm = 1.9212 # pbe0
e0 = -911.5521164109 * Ry2eV
e_n1 = -911.6816801886 * Ry2eV + (-0.06010319) * Ry2eV * 2
e_p1 = -911.4850239942 * Ry2eV + (0.10517428) * Ry2eV * 2
e_n2 =  -911.7312241524 * Ry2eV + (-0.07458343) * Ry2eV * 2
e_p2 = -911.3136421332 * Ry2eV + (0.25555818) * Ry2eV * 2


# 6x6 CBVN in-plane jdftx
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781  - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 # pbe0
#cbm = 1.9212 # pbe0
#e0 = -454.6127977032838317 * Ry2eV * 2
#e_n1 = -454.8309898573902501 * Ry2eV * 2 + (0.04749697) * Ry2eV * 2
#e_p1 = -454.5047200242801182 * Ry2eV * 2 + (0.04226736) * Ry2eV * 2
#e_n2 =  * Ry2eV * 2 + (-0.07458343) * Ry2eV * 2
#e_p2 =  * Ry2eV * 2 + (0.25555818) * Ry2eV * 2



# 6x6 CBVN out-of-plane qe+jdftx
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781  - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 # pbe0
#cbm = 1.9212 # pbe0
#e0 = -911.6885743725 * Ry2eV
#e_n1 = -911.6822001847 * Ry2eV + (-0.06000355) * Ry2eV * 2
#e_p1 = -911.5420129342 * Ry2eV + (0.10546195) * Ry2eV * 2
#e_n2 =  -911.7317280418 * Ry2eV + (-0.07420310) * Ry2eV * 2
#e_p2 = -911.3137821287 * Ry2eV + (0.25593845) * Ry2eV * 2



# 6x6 CB-hBN
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781 - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 # pbe0
#cbm = 1.9212 # pbe0
#e0 = -931.9875450170 * Ry2eV
#e_n1 = -911.6822001847 * Ry2eV + (-0.06000355) * Ry2eV * 2
#e_p1 = -931.9356209946 * Ry2eV + (0.10704068) * Ry2eV * 2
#e_n2 =  -911.7317280418 * Ry2eV + (-0.07420310) * Ry2eV * 2
#e_p2 = -911.3137821287 * Ry2eV + (0.25593845) * Ry2eV * 2


# 6x6 NBVN qe+jdftx
#vbm = -3.4945 - 0.171369164 * Ry2eV #pbe
#cbm = 1.1781 - 0.171369164 * Ry2eV #pbe
#vbm = -5.4414 # pbe0
#cbm = 1.9212 # pbe0
#e0 = -920.1466536071 * Ry2eV
#e_n1 = -920.1719185150 * Ry2eV + (-0.06031619) * Ry2eV * 2
#e_p1 = -919.9415051256 * Ry2eV + (0.10601984) * Ry2eV * 2
#e_n2 = -920.1082435283 * Ry2eV + (-0.08439696) * Ry2eV * 2
#e_p2 = -919.6971911899 * Ry2eV + (0.25703454) * Ry2eV * 2


plt.axvline(e0-e_p1-vbm, color='k', ls='dotted')
plt.axvline(e_n1-e0-vbm, color='k', ls='dotted')


x = np.linspace(0, cbm-vbm)

e0 = e0 + x*0
e_n1 = neg(e_n1, x, vbm, 1)
e_p1 = pos(e_p1, x, vbm, 1)
#e_n2 = neg(e_n2, x, vbm, 2)
#e_p2 = pos(e_p2, x, vbm, 2)


plt.plot(x, e0, label="q=0")
plt.plot(x, e_n1, label="q=-1")
plt.plot(x, e_p1, label="q=+1")
#plt.plot(x, e_n2, label="q=-2")
#plt.plot(x, e_p2, label="q=+2")
plt.legend()
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.xlabel("$\epsilon_{fermi}-VBM$ (eV)")
plt.ylabel("Formation Energy (eV)")
plt.show()