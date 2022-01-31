#!/usr/local/bin/python

from pylab import *
import argparse
import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_fwf("trial.dat", header=None)
#data2 = pd.read_fwf("problem1b2.dat", header=None)
#data3 = pd.read_fwf("problem1b3.dat", header=None)


plt.figure(figsize=(10,10), dpi=120)
plt.plot(data1[1], data1[3], color='b', label='(1,0)')
#plt.plot(data2[2], data2[3], color='r', label='(2,0)')
#plt.plot(data3[2], data3[3], color='g', label='(0,3)')
#plt.xlim(-5, 10)
#plt.ylim(-5, 10)
plt.title("Problem 1b")
plt.xlabel("x")
plt.ylabel("y")
plt.legend(loc="best")
plt.show()





