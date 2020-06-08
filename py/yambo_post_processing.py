#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
from constants import *
from configuration_for_plot import config_plot


class yambo_post_processing:
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
                index_k_cbm = np.where(self.eigenE[0:self.nk, vbm] + \
                                        self.corr[0:self.nk, vbm] == \
                                    np.amin(self.eigenE[0:self.nk, vbm] + \
                                            self.corr[0:self.nk, vbm])
                                    )[0][0]
                index_k_vbm = np.where(self.eigenE[0:self.nk, vbm-1] + \
                                    self.corr[0:self.nk, vbm-1] == \
                                np.amax(self.eigenE[0:self.nk, vbm-1] + \
                                    self.corr[0:self.nk, vbm-1])
                                    )[0][0]
            else:
                index_k_cbm = np.where(self.eigenE[self.nk:self.nk*2, vbm] + \
                                    self.corr[self.nk:self.nk*2, vbm] == \
                                np.amin(self.eigenE[self.nk:self.nk*2, vbm] + \
                                    self.corr[self.nk:self.nk*2, vbm])
                                    )[0][0]
                index_k_vbm = np.where(self.eigenE[self.nk:self.nk*2, vbm-1] + \
                                    self.corr[self.nk:self.nk*2, vbm-1] == \
                            np.amax(self.eigenE[self.nk:self.nk*2, vbm-1] + \
                                    self.corr[self.nk:self.nk*2, vbm-1])
                                    )[0][0]
        else:
            self.eigenE = self.data[:, 2].reshape((self.nk, self.nbnd))
            self.corr = self.data[:, 3].reshape((self.nk, self.nbnd))
            self.direct_gap = self.eigenE[:, vbm] - self.eigenE[:, vbm-1] + \
                                self.corr[:, vbm] - self.corr[:, vbm-1]
            self.indirect_gap = np.amin(self.eigenE[:, vbm] + \
                                    self.corr[:, vbm]) - \
                                np.amax(self.eigenE[:, vbm-1] + \
                                    self.corr[:, vbm-1])
            index_k_cbm = np.where(self.eigenE[:, vbm] + \
                                    self.corr[:, vbm] == \
                            np.amin(self.eigenE[:, vbm] + \
                                    self.corr[:, vbm])
                                    )[0][0]
            index_k_vbm = np.where(self.eigenE[:, vbm-1] + \
                                    self.corr[:, vbm-1] == \
                            np.amax(self.eigenE[:, vbm-1] + \
                                    self.corr[:, vbm-1])
                                    )[0][0]
        
        sys.stdout.write("CBM is at No.{} K point; CBM = {} eV\n"\
                        .format(index_k_cbm, \
                        self.eigenE[index_k_cbm, vbm] + \
                        self.corr[index_k_cbm, vbm]))
        sys.stdout.write("VBM is at No.{} K point; VBM = {} eV\n"\
                        .format(index_k_vbm, \
                        self.eigenE[index_k_vbm, vbm-1] + \
                        self.corr[index_k_vbm, vbm-1]))
        sys.stdout.write("Bandgap = {} eV\n".format(self.indirect_gap))
        sys.stdout.flush()
        return(self.corr)

if __name__ == "__main__":

    config_plot()
    path = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs"
    bands = ["nk18_nbnd100_qe6.1_yambo4.4", "nk18_nbnd200_qe6.1_yambo4.4",
    "nk18_nbnd300_qe6.1_yambo4.4", "nk18_nbnd400_qe6.1_yambo4.4"]
    rys = ["data-1ry/o-all_Bz.qp", "data-3ry/o-all_Bz.qp", "data-5ry/o-all_Bz.qp", "data-8ry/o-all_Bz.qp", "data-12ry/o-all_Bz.qp", "data-15ry/o-all_Bz.qp", "data-20ry/o-all_Bz.qp"]

    gap = np.zeros((len(rys),len(bands)), dtype=float)
    label1 = ["1 Ry", "3 Ry", "5 Ry", "8 Ry", "12 Ry", "15ry", "20ry"]
    label2 = ["BndsRnXp=100", "BndsRnXp=200", "BndsRnXp=300", "BndsRnXp=400"]
    for i, ry in enumerate(rys):
        for j, band in enumerate(bands):
            dir = os.path.join(path, band, ry)
            yambo = yambo_post_processing(dir)
            correction = yambo.read_corr(vbm=4)
            gap[i, j] = yambo.indirect_gap

    for i in range(len(bands)):
        plt.plot([1, 3, 5, 8, 12, 15, 20], gap[:, i], label=label2[i])
    plt.scatter([1, 3, 5, 8, 12, 15, 20], gap[:, 0])
    plt.scatter([1, 3, 5, 8, 12, 15, 20], gap[:, 1])
    plt.scatter([1, 3, 5, 8, 12, 15, 20], gap[:, 2])
    plt.scatter([1, 3, 5, 8, 12, 15, 20], gap[:, 3])
    plt.legend()
    plt.xlabel("NGsBlkXp (Ry)")
    plt.ylabel("$E_{gap (K → Γ)}$ (eV)")
    plt.title("GW@PBE0 (⍺=0.41)")
    plt.show()


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


    """
    dir = "/home/likejun/work/pristine_hbn/1x1/exx0.41/scf-gs/nk18_nbnd300_qe6.1_yambo4.4/data-3ry/o-all_Bz.qp"
    yambo = yambo_post_processing(dir)
    yambo.read_corr(4)
    """
