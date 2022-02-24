import matplotlib.pyplot as plt
import numpy as np
import math

M=10**11*4.30091*10**(-3)*(1.02201*10**(-6))**2
# An "interface" to matplotlib.axes.Axes.hist() method
filename = "init.dat"
x = np.loadtxt(filename, usecols = (1), skiprows=1)
x2 = np.linspace(0, 20000, 100000)
rho = 3.0/(4.0*math.pi)*1500**(-3)*(1+(x2/1500)**2)**(-5.0/2.0)*4*math.pi*x2**2 # since mass density, treat total mass as 1
fig1, ax1 = plt.subplots()
n, bins, patches = ax1.hist(x, bins='auto', density=True, color='#0504aa', alpha=0.75)
ax1.plot(x2, rho, 'r', label = 'theoretical distribution of mass density')
ax1.legend()
ax1.set_xlabel('radial distance from the origin, kparsec')
ax1.set_ylabel('Mass density (unit = 10^11 solarmass)')
ax1.set_title('mass distribution of the radial distance')
ax1.set_xlim(0,50000)

ax1.grid(True)
fig1.savefig("mass_distr.png")

