import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use("/home/likejun/work/github/plot_tools/styles/wamum")
S = [ 180, 180, 192 ]            #Sample count (fftbox from log file)
L = [ 1, 1, 30 ]    #Lattice vector lengths (bohrs)
#S = [ 140, 140, 96 ] # redo_wufeng_example
#L = [ 1, 1, 20  ]    #Lattice vector lengths (bohrs)


z = np.arange(S[2])*L[2]/S[2]  #z-coordinates of grid points

def readPlanarlyAveraged(fileName, S):
    out = np.fromfile(fileName, dtype=np.float64)
    out = np.reshape(out, S)   #Reshape data
    out = np.mean(out, axis=1) #y average
    out = np.mean(out, axis=0) #x average
    return out

directory = "/home/likejun/ctl/cb/cb/redo_tyler_example/repeat_2-JDFTx_wo_Cut/q+1"
dir_dtot_charge_p1 = os.path.join(directory, "charge.d_tot")
dir_dtot_neutral = os.path.join(directory, "neutral.d_tot")
dir_prist = os.path.join(directory, "prist.d_tot")
dtot_charge_p1 = readPlanarlyAveraged(dir_dtot_charge_p1, S)
dtot_neutral = readPlanarlyAveraged(dir_dtot_neutral, S)
prist_Dtot = readPlanarlyAveraged(dir_prist, S)
#dShift = dir_dtot_charge_p1 - dtot_neutral

#print("VBMshift =", dShift[0]*27.2114)  #Report VBM shift in eV
#plt.plot(z*0.529, dShift*27.2114, label="Dvac - Dtot");     #Plot with unit conversions
plt.plot(z*0.529, dtot_charge_p1*27.2114, label="Dtot (q=+1)")
plt.plot(z*0.529, dtot_neutral*27.2114, label="Dtot (q=0)")
plt.plot(z*0.529, prist_Dtot*27.2114, label="Dtot (prist)")
plt.xlabel('z [Angstroms]')
plt.ylabel('Potential [eV]')
plt.xlim([0, L[2]*0.529])            #Select first half (top surface only)
plt.legend()
plt.show()
