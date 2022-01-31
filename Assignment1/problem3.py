#!/usr/local/bin/python

from pylab import *
import argparse
import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_fwf("Marsanalytic.dat", header=None)
data2 = pd.read_fwf("planetorbitdat/Mars.dat", header=None)


plt.figure(figsize=(10,10), dpi=180, linewidth=0.8)
plt.rcParams.update({'font.size': 6})
plt.plot(data1[0], data1[1], color='palegreen', linewidth=0.5, label='Analytic orbit')
plt.plot(data2[1], data2[3], color='dodgerblue', linewidth=0.5, label='Leapfrog simulation')
ax = plt.gca()
ax.set_facecolor('black')
plt.title("Problem 3b for Mars's orbit")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc="best")
plt.savefig("Problem3Mars.png")




