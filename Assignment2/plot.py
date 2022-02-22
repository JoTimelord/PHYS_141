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
ax1.set_ylabel('Mass density (unit = 10^11 solarmass)')
ax1.set_title('(normalized) mass density distribution of the radial distance')
ax1.grid(True)
fig1.savefig("mass_distr.png")

E = np.loadtxt(filename, usecols = (9), skiprows=1)
E = E/1000000 # put the energy at (Mm/s)^2
E2 = np.linspace(0, -0.3, 100000)
f = 24*math.sqrt(2)/(7*math.pi**3)*(4.302*10**5)**(-4)*1.5**2*(-E2)**(3.5)*np.sqrt(-2*E2)
#f = 24*math.sqrt(2)/(7*math.pi**3)*(-E2)**(3.5)*np.sqrt(2*-E2)
fig2, ax2 = plt.subplots()
n, bins, patches = ax2.hist(E, density=False, bins='auto', color='red', alpha=0.75)
ax2.plot(E2, f, 'b', label = 'theoretical distribution of energy')
ax2.legend()
ax2.set_xlabel('energy distribution of points, in 10^6(km/s)^2')
ax2.set_ylabel('Frequency')
ax2.set_title('energy distribution of points')
ax2.grid(True)
fig2.savefig("energy_distr.png")


