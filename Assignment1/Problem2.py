#!/usr/local/bin/python

from pylab import *
import argparse
import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_fwf("Jupiter.dat", header=None)
data2 = pd.read_fwf("Neptune.dat", header=None)
data3 = pd.read_fwf("Mars.dat", header=None)
data4 = pd.read_fwf("Earth.dat", header=None)
data5 = pd.read_fwf("Pluto.dat", header=None)
data6 = pd.read_fwf("Venus.dat", header=None)
data7 = pd.read_fwf("Saturn.dat", header=None)
data8 = pd.read_fwf("Uranus.dat", header=None)
data9 = pd.read_fwf("Mercury.dat", header=None)

plt.figure(figsize=(20,20), linewidth=0.9, dpi=180)
plt.rcParams.update({'font.size': 6})
plt.plot(data1[1], data1[3], color='slateblue', linewidth=0.5, label='Jupiter')
plt.plot(data2[1], data2[3], color='gold', linewidth=0.5, label='Neptune')
plt.plot(data3[1], data3[3], color='aqua', linewidth=0.5, label='Mars')
plt.plot(data4[1], data4[3], color='powderblue', linewidth=0.5, label='Earth')
plt.plot(data5[1], data5[3], color='darkgreen', linewidth=0.5, label='Pluto')
plt.plot(data6[1], data6[3], color='sienna', linewidth=0.5, label='Venus')
plt.plot(data7[1], data7[3], color='steelblue', linewidth=0.5, label='Saturn')
plt.plot(data8[1], data8[3], color='chocolate', linewidth=0.5, label='Uranus')
plt.plot(data9[1], data9[3], color='olivedrab', linewidth=0.5, label='Mercury')
ax = plt.gca()
ax.set_facecolor('black')
plt.title("Problem 2: an ugly composite of all planetary orbits")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc="best")
plt.savefig("planetorbitplot.png")


