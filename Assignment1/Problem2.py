#!/usr/local/bin/python

from pylab import *
import argparse
import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_fwf("Jupiter.dat", header=None)

print(data1[2][1], data1[4][1])
plt.figure(figsize=(10,10), dpi=120)
plt.plot(data1[2]/1000000, data1[4]/1000000, color='b', label='Jupiter')
#plt.xlim(-300, 200)
#plt.ylim(-300, 200)
plt.title("Problem 1b")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc="best")
plt.show()




