import matplotlib.pyplot as plt
import numpy as np
import math

# An "interface" to matplotlib.axes.Axes.hist() method
filename = "backup.dat"

E = np.loadtxt(filename, usecols = (9), skiprows=1)
E = E # put the energy at (Mm/s)^2
E2 = np.linspace(0, -1, 100000)
f = 24*math.sqrt(2)/(7*math.pi**3)*(-E2)**(3.5)*np.sqrt(2*-E2)
fig2, ax2 = plt.subplots()
n, bins, patches = ax2.hist(E, density=True, bins='auto', color='red', alpha=0.75)
ax2.plot(E2, f, 'b', label = 'theoretical distribution of energy')
ax2.legend()
ax2.set_xlabel('energy distribution of points')
ax2.set_ylabel('Frequency')
ax2.set_title('energy distribution of points')
ax2.grid(True)
fig2.savefig("energy_distr_back_up.png")


