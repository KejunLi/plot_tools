#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
from constants import *
from configuration_for_plot import config_plot


class yambo_post_processing:
    """
    ============================================================================
    +   1. Constructor
    +   Attributes:
    +   self.fname (the specific directory to the file)
    +   self.dir (the directory which contains the file)
    +   self.yambo_out (read the file)
    +   self.lines (lines in the file)
    +   self.data (data in the file)
    +   self.spinpol (spin polarization condition)
    +   self.nk (number of k points)
    +   self.kpoints (indies of k points included)
    +   self.nbnd (number of bands (Kohn-Sham states))
    +
    +   No return
    ============================================================================
    +   2. Method read_corr(self, index_of_vbm)
    +   Attributes:
    +   self.eigenE (eigenenergies, eV)
    +   self.corr (quasi-particle corrections)
    +   self.direct_gap (direct gap at each k point)
    +   self.indirect_gap (indirect gap between CBM and VBM)
    +   self.cbm (value of CBM)
    +   self.vbm (value of VBM)
    +
    +   No return
    ============================================================================
    """
    def __init__(self, dir_f):
        """
        initial method or constructor for initialization
        """
        self.fname = dir_f
        self.dir = os.path.dirname(dir_f)
        if self.dir == "":
            self.dir = "."

        try:
            self.yambo_out = open(self.fname, "r")
        except:
            raise IOError("Fail to open {}".format(self.fname))

        self.lines = self.yambo_out.readlines()
        self.data = np.genfromtxt(dir_f, dtype=None)
        self.spinpol = False

        for line in self.lines:
            if "QP @ K" in line:
                k_i = int(re.findall(r"[+-]?\d+", line)[0]) # initial k point
                k_f = int(re.findall(r"[+-]?\d+", line)[1]) # final k point
                band_i = int(re.findall(r"[+-]?\d+", line)[2]) # initial band
                band_f = int(re.findall(r"[+-]?\d+", line)[3]) # final band
                self.nk = k_f - k_i + 1
                self.kpoints = np.array(range(self.nk))
                self.nbnd = band_f - band_i + 1
            if "Spin_Pol" in line:
                self.spinpol = True
        sys.stdout.write("\rYambo\n")
        sys.stdout.write("Quasi-particle correction range:\n")
        sys.stdout.write("K kpoints: {}\n".format(self.kpoints))
        sys.stdout.write("Bands: {}\n".format((band_i, band_f)))
        sys.stdout.flush()

    def read_corr(self, vbm):
        if self.spinpol:
            eigenE = self.data[:, 2].reshape((self.nk, self.nbnd*2))
            corr = self.data[:, 3].reshape((self.nk, self.nbnd*2))

            self.eigenE = np.zeros((self.nk*2, self.nbnd), dtype=float)
            self.corr = np.zeros((self.nk*2, self.nbnd), dtype=float)
            eigenE_up = np.zeros((self.nk, self.nbnd), dtype=float)
            eigenE_dn = np.zeros((self.nk, self.nbnd), dtype=float)
            corr_up = np.zeros((self.nk, self.nbnd), dtype=float)
            corr_dn = np.zeros((self.nk, self.nbnd), dtype=float)
            self.direct_gap = np.zeros(self.nk*2, dtype=float)

            for i in range(self.nbnd*2):
                if i % 2 == 0: # spin up values
                    eigenE_up[:, i//2] = eigenE[:, i]
                    corr_up[:, i//2] = corr[:, i]
                else: # spin down values
                    eigenE_dn[:, i//2] = eigenE[:, i]
                    corr_dn[:, i//2] = corr[:, i]

            self.eigenE[0:self.nk, :] = eigenE_up
            self.eigenE[self.nk:self.nk*2, :] = eigenE_dn
            self.corr[0:self.nk, :] = corr_up
            self.corr[self.nk:self.nk*2, :] = corr_dn

            self.direct_gap = self.eigenE[:, vbm] - self.eigenE[:, vbm-1] + \
                                self.corr[:, vbm] - self.corr[:, vbm-1]

            indirect_gap_up = np.amin(self.eigenE[0:self.nk, vbm] + \
                                    self.corr[0:self.nk, vbm]) - \
                            np.amax(self.eigenE[0:self.nk, vbm-1] + \
                                    self.corr[0:self.nk, vbm-1])
            indirect_gap_dn = np.amin(self.eigenE[self.nk:self.nk*2, vbm] + \
                                self.corr[self.nk:self.nk*2, vbm]) - \
                            np.amax(self.eigenE[self.nk:self.nk*2, vbm-1] + \
                                self.corr[self.nk:self.nk*2, vbm-1])
            if indirect_gap_up > 0 and indirect_gap_dn > 0:
                self.indirect_gap = min(indirect_gap_up, indirect_gap_dn)
            elif indirect_gap_up < 0 and indirect_gap_dn < 0:
                raise ValueError("Negative gap")
            else:
                self.indirect_gap = max(indirect_gap_up, indirect_gap_dn)

            # find the index of k points where vbm and cbm are
            if self.indirect_gap == indirect_gap_up:
                self.cbm = np.amin(self.eigenE[0:self.nk, vbm] + \
                                    self.corr[0:self.nk, vbm])
                self.vbm = np.amax(self.eigenE[0:self.nk, vbm-1] + \
                                    self.corr[0:self.nk, vbm-1])
                index_k_cbm = np.where(self.eigenE[0:self.nk, vbm] + \
                                    self.corr[0:self.nk, vbm] == self.cbm
                                    )[0][0]
                index_k_vbm = np.where(self.eigenE[0:self.nk, vbm-1] + \
                                    self.corr[0:self.nk, vbm-1] == self.vbm
                                    )[0][0]
            else:
                self.cbm = np.amin(self.eigenE[self.nk:self.nk*2, vbm] + \
                                    self.corr[self.nk:self.nk*2, vbm])
                self.vbm = np.amax(self.eigenE[self.nk:self.nk*2, vbm-1] + \
                                    self.corr[self.nk:self.nk*2, vbm-1])
                index_k_cbm = np.where(self.eigenE[self.nk:self.nk*2, vbm] + \
                                    self.corr[self.nk:self.nk*2, vbm] == \
                                    self.cbm
                                    )[0][0]+self.nk
                index_k_vbm = np.where(self.eigenE[self.nk:self.nk*2, vbm-1] + \
                                    self.corr[self.nk:self.nk*2, vbm-1] == \
                                    self.vbm
                                    )[0][0]+self.nk
        else:
            self.eigenE = self.data[:, 2].reshape((self.nk, self.nbnd))
            self.corr = self.data[:, 3].reshape((self.nk, self.nbnd))
            self.direct_gap = self.eigenE[:, vbm] - self.eigenE[:, vbm-1] + \
                                self.corr[:, vbm] - self.corr[:, vbm-1]
            self.cbm = np.amin(self.eigenE[:, vbm] + self.corr[:, vbm])
            self.vbm = np.amax(self.eigenE[:, vbm-1] + self.corr[:, vbm-1])
            self.indirect_gap = self.cbm - self.vbm
            index_k_cbm = np.where(self.eigenE[:, vbm] + \
                                    self.corr[:, vbm] == self.cbm
                                    )[0][0]
            index_k_vbm = np.where(self.eigenE[:, vbm-1] + \
                                    self.corr[:, vbm-1] == self.vbm
                                    )[0][0]

        sys.stdout.write("CBM is at No.{} K point; CBM = {} eV\n"\
                        .format(index_k_cbm, self.cbm))
        sys.stdout.write("VBM is at No.{} K point; VBM = {} eV\n"\
                        .format(index_k_vbm, self.vbm))
        sys.stdout.write("Bandgap = {} eV\n".format(self.indirect_gap))
        sys.stdout.flush()
        return(self.corr)

if __name__ == "__main__":

    config_plot()
#===============================================================================
    path = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs"
    bands = ["nk18_nbnd100_qe6.1_yambo4.4", "nk18_nbnd200_qe6.1_yambo4.4",
            "nk18_nbnd300_qe6.1_yambo4.4", "nk18_nbnd400_qe6.1_yambo4.4",
            "nk18_nbnd500_qe6.1_yambo4.4", "nk18_nbnd600_qe6.1_yambo4.4",
            "nk18_nbnd700_qe6.1_yambo4.4"]
    rys = ["data-1ry/o-all_Bz.qp", "data-3ry/o-all_Bz.qp",
            "data-5ry/o-all_Bz.qp", "data-8ry/o-all_Bz.qp",
            "data-12ry/o-all_Bz.qp", "data-15ry/o-all_Bz.qp",
            "data-20ry/o-all_Bz.qp", "data-25ry/o-all_Bz.qp", "data-30ry/o-all_Bz.qp"]
    x = [1, 3, 5, 8, 12, 15, 20, 25, 30]
    label1 = ["1 Ry", "3 Ry", "5 Ry", "8 Ry", "12 Ry", "15ry", "20ry", "25ry", "30Ry"]
    label2 = ["BndsRnXp=100", "BndsRnXp=200", "BndsRnXp=300", "BndsRnXp=400",
            "BndsRnXp=500", "BndsRnXp=600", "BndsRnXp=700"]

    """
    path = "/home/likejun/work/pristine_hbn/1x1/qe_convergence/pbe/bands"
    bands = ["dnbnd_800", "dnbnd_900",]
    rys = ["data-1ry/o-all_Bz.qp", "data-3ry/o-all_Bz.qp",
            "data-5ry/o-all_Bz.qp", "data-8ry/o-all_Bz.qp",
            "data-12ry/o-all_Bz.qp", "data-15ry/o-all_Bz.qp",
            "data-20ry/o-all_Bz.qp"]#, "data-25ry/o-all_Bz.qp"]
    x = [1, 3, 5, 8, 12, 15, 20]
    label1 = ["1 Ry", "3 Ry", "5 Ry", "8 Ry", "12 Ry", "15ry", "20ry"]#, "25ry"]
    label2 = ["BndsRnXp=800", "BndsRnXp=900"]
    """
    gap = np.zeros((len(rys),len(bands)), dtype=float)
    cbm = np.zeros((len(rys),len(bands)), dtype=float)
    vbm = np.zeros((len(rys),len(bands)), dtype=float)
    Ef = -5.441420 # fermi level of exx=0.41
    #Ef = -4.658698 # fermi level of exx=0.25
    Vvac = 2.3355 # electrostatic potential
    for i, ry in enumerate(rys):
        for j, band in enumerate(bands):
            dir = os.path.join(path, band, ry)
            print(dir)
            yambo = yambo_post_processing(dir)
            correction = yambo.read_corr(vbm=4)
            gap[i, j] = yambo.indirect_gap
            cbm[i, j] = yambo.cbm - Vvac + Ef
            vbm[i, j] = yambo.vbm - Vvac + Ef
    print(gap)
    for i in range(len(bands)):
        plt.plot(x, gap[:, i], label=label2[i])
        plt.scatter(x, gap[:, i])

    #for i in range(len(bands)):
    #    plt.plot(x, cbm[:, i], label=label2[i])
    #    plt.scatter(x, cbm[:, i])
    left, right = plt.xlim()
    #plt.hlines(0.0, left, right, linestyles="dashed", linewidth=0.8)
    plt.legend()
    plt.xlabel("NGsBlkXp (Ry)")
    #plt.ylabel("$VBM (K)$ (eV)")
    #plt.ylabel("$CBM (Γ)$ (eV)")
    plt.ylabel("$E_{gap (K → Γ)}$ (eV)")
    #plt.title("GW@PBE\n nk=18x18x1, ecutwfc=60 Ry") # PBE
    plt.title("GW@PBE0 (⍺=0.41) \n nk=18x18x1, ecutwfc=60 Ry, nqx=6x6x1") # exx=0.41
    #plt.title("GW@PBE0 \n nk=18x18x1, ecutwfc=60 Ry, nqx=6x6x1") # exx=0.25

    plt.xlim()
    plt.show()
#===============================================================================

    """
    config_plot()
    path = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd300_qe6.1_yambo4.4"
    rys = ["data-1ry/o-all_Bz.qp", "data-3ry/o-all_Bz.qp", "data-5ry/o-all_Bz.qp", "data-8ry/o-all_Bz.qp", "data-12ry/o-all_Bz.qp", "data-15ry/o-all_Bz.qp", "data/o-all_Bz.qp"]

    gap = np.zeros(7, dtype=float)
    label1 = ["1 Ry", "3 Ry", "5 Ry", "8 Ry", "12 Ry", "15ry", "20ry"]
    for i, ry in enumerate(rys):
        dir = os.path.join(path, ry)
        yambo = yambo_post_processing(dir)
        correction = yambo.read_corr(vbm=4)
        gap[i] = yambo.indirect_gap

    print(gap)

    plt.plot([1, 3, 5, 8, 12, 15, 20], gap)
    plt.xlabel("NGsBlkXp (Ry)")
    plt.ylabel("$E_{gap (K → Γ)}$ (eV)")
    plt.title("GW@PBE0 (⍺=0.41)")
    plt.xlim(-1.7,21)
    plt.show()
    """
#===============================================================================

    """
    dir = "/home/likejun/data/o-all_Bz.qp"
    yambo = yambo_post_processing(dir)
    yambo.read_corr(4)
    print(yambo.eigenE[0,4], yambo.corr[0,4])
    print(yambo.eigenE[-1,3], yambo.corr[-1,3])
    """
#===============================================================================
