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
        component_1 = np.exp(k*r*1j)
        component_2 = (1-k*r*1j)*(3*costheta**2-1)/r**3
        component_3 = (k*sintheta)**2/r
        scatt = component_1 * (component_2 + component_3)
        Scattering += scatt
        #print(Scattering)
    return(np.real(Scattering), np.imag(Scattering))
    #return(Scattering)

def gen_dipole_array(dim):
    """
    function used for generating an array centered at [0, 0, 0]
    """
    number_of_dipoles = dim**2
    val = np.array([-dim, -dim, 0])
    dipole_positions = np.full((number_of_dipoles, 3), (val+val%2)/2)
    for i in range(dim):
        for j in range(dim):
            if dim%2 == 0:
                dipole_positions[i*dim+j, 0] += i + 0.5
                dipole_positions[i*dim+j, 1] += j + 0.5
            else:
                dipole_positions[i*dim+j, 0] += i
                dipole_positions[i*dim+j, 1] += j
    return(dipole_positions[~(dipole_positions==0).all(1)])
    

origin = np.array([0, 0, 0])
dim = 100 # array dimension
dips = gen_dipole_array(dim)/10**9 * 300 # dipole positions
#print(dips)

dip_orientation = np.array([1, 0, 0])
k = 2*np.pi/(np.arange(450, 470, 0.1)*10**(-9))
Re_S = np.zeros(len(k))
Im_S = np.zeros(len(k))
for i in range(len(k)):
    Re_S[i] = scattering(origin, dips, dip_orientation, k[i])[0]
    Im_S[i] = scattering(origin, dips, dip_orientation, k[i])[1]

plt.rcParams['figure.figsize'] = (12, 9)
plt.rcParams['figure.dpi'] = 120
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.plot(k, Re_S, color="tab:red", label="Re(S)")
plt.plot(k, Im_S, color="tab:blue", label="Im(S)")
plt.xlabel("k")
plt.ylabel("S")
plt.legend()

plt.show()