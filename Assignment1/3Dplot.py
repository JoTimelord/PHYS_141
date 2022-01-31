#!/usr/local/bin/python

from pylab import *
import argparse
from mpl_toolkits import mplot3d
import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_fwf("problem1a1.dat", header=None)
data2 = pd.read_fwf("problem1a2.dat", header=None)
data3 = pd.read_fwf("problem1a3.dat", header=None)

fig = plt.figure(figsize=(10,10), dpi=180)
ax = plt.axes(projection ='3d')
ax.set_title('3D line plot geeks for geeks')
ax.plot3D(data1[1], data1[3], data1[0],'green', label='(1,0)')
ax.plot3D(data2[1], data2[3], data2[0],'blue', label='(2,0)')
ax.plot3D(data3[1], data3[3], data3[0],'red', label='(0,3)')
plt.title("3D plot for linear")
plt.xlabel("x")
plt.ylabel("v")
ax.set_zlabel("t")
plt.legend(loc="best")
plt.savefig('3Dplot1')






