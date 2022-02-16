import matplotlib.pyplot as plt
import numpy as np
import math

# An "interface" to matplotlib.axes.Axes.hist() method
filename = "1a.dat"
x = np.loadtxt(filename, usecols = (1), skiprows=1)
x2 = np.linspace(0, 17, 100000)
rho = 3.0/(4.0*math.pi)*1.5**(-3)*(1+(x2/1.5)**2)**(-5.0/2.0)*4*math.pi*x2**2 # since mass density, treat total mass as 1
fig1, ax1 = plt.subplots()
n, bins, patches = ax1.hist(x, bins='auto', density=True, color='#0504aa', alpha=0.75)
ax1.plot(x2, rho, 'r', label = 'theoretical distribution of mass density')
ax1.legend()
ax1.set_xlabel('radial distance from the origin, kparsec')
ax1.set_ylabel('Frequency')
ax1.set_title('mass density distribution of the radial distance')
ax1.grid(True)
ax1.set_xlim(0, 15.0)
ax1.set_ylim(0, 1.0)
fig1.savefig("mass_distr.png")

E = np.loadtxt(filename, usecols = (9), skiprows=1)
fig2, ax2 = plt.subplots()
n, bins, patches = ax2.hist(E, bins='auto', density=True, color='red', alpha=0.75)
# ax2.plot(x2, rho, 'r', label = 'theoretical distribution of mass density')
# ax2.legend()
ax2.set_xlabel('energy distribution of points, in (km/s)^2')
ax2.set_ylabel('Frequency')
ax2.set_title('energy distribution of points')
ax2.grid(True)
ax2.set_xlim(0, 15.0)
ax2.set_ylim(0, 1.0)
fig2.savefig("energy_distr.png")


