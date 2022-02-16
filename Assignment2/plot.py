import matplotlib.pyplot as plt
import numpy as np

# An "interface" to matplotlib.axes.Axes.hist() method
filename = "1a.dat"
x = np.loadtxt(filename, usecols = (1), skiprows=1)
n, bins, patches = plt.hist(x, bins='auto', density=True, color='#0504aa', alpha=0.75)
plt.xlabel('radial distance from the origin, kparsec')
plt.ylabel('Frequency')
plt.title('mass density distribution of the radial distance')
plt.grid(True)
plt.xlim(0, 15.0)
plt.ylim(0, 1.0)
plt.savefig("mass_distr.png")
