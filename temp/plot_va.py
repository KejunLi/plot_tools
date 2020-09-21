import numpy as np
import matplotlib.pyplot as plt

def scattering(origin, dipole_positions, dipole_orientation, k):
    Scattering = 0
    for i, dip in enumerate(dipole_positions):
        r_vec = origin - dip # displacement
        r = np.sqrt(np.dot(r_vec, r_vec)) # distance
        e = r_vec/r # direction vector
        costheta = np.dot(e, dipole_orientation)
        sintheta = np.sqrt(1 - costheta**2)
        scatt = np.exp(np.imag(k*r*1j)) * \
            (((1-k*r*1j)*(3*costheta**2-1))/r**3+(k*sintheta)**2/r)
        Scattering += scatt
        #print(Scattering)
    return(np.real(Scattering), np.imag(Scattering))

def gen_dip(dimension):
    number_of_dipoles = dimension**2
    dipole_positions = np.zeros((number_of_dipoles, 3))
    for i in range(dimension):
        for j in range(dimension):
            dipole_positions[i*dimension+j, 0] = i + 1
            dipole_positions[i*dimension+j, 1] = j + 1
    return(dipole_positions)
    

origin = np.array([0, 0, 1])
dim = 10
dips = gen_dip(dim)
dip_orientation = [0, 1, 0]
k = np.arange(1, 5.0, 0.1)
Re_S = np.zeros(len(k))
Im_S = np.zeros(len(k))
for i in range(len(k)):
    Re_S[i] = scattering(origin, dips, dip_orientation, k[i])[0]
    Im_S[i] = scattering(origin, dips, dip_orientation, k[i])[1]


plt.plot(k, Re_S)
plt.plot(k, Im_S)
plt.show()

