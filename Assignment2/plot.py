import matplotlib.pyplot as plt
import numpy as np

# An "interface" to matplotlib.axes.Axes.hist() method
filename = "1a.dat"
x = np.loadtxt(filename, usecols = (1), skiprows=1)
n, bins, patches = plt.hist(x, bins='auto', density=True, color='#0504aa', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('radial distance from the origin, kparsec')
plt.ylabel('Frequency')
plt.title('distribution of the radial distance')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.savefig("mass_distr.png")
