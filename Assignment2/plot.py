import matplotlib.pyplot as plt
import numpy as np
import math

# An "interface" to matplotlib.axes.Axes.hist() method
filename = "1a.dat"
x = np.loadtxt(filename, usecols = (1), skiprows=1)
x2 = np.linspace(0, 17, 100000)
rho = 3.0/(4.0*math.pi)*1.5**(-3)*(1+(x2/1.5)**2)**(-5.0/2.0)*4*math.pi*x2**2 # since mass density, treat total mass as 1
n, bins, patches = plt.hist(x, bins='auto', density=True, color='#0504aa', alpha=0.75)
plt.plot(x2, rho, 'r', label = 'theoretical distribution of mass density')
plt.legend()
plt.xlabel('radial distance from the origin, kparsec')
plt.ylabel('Frequency')
plt.title('mass density distribution of the radial distance')
plt.grid(True)
plt.xlim(0, 15.0)
plt.ylim(0, 1.0)
plt.savefig("mass_distr.png")
