#!/usr/local/bin/python
#***********************************************************************#
#*                      TRAJECTORY ANIMATION                           *#
#*                                                                     *#
#* This program reads the X-Y trajectory of a particle an generates    *#
#* an animation of the motion.                                         *#
#*                                                                     *#
#***********************************************************************#
 # Adapted from AUTHOR: FELIPE GONZALEZ CATALDO, September 2018.
 # Referenced: https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30
from pylab import *
import imageio
import os

fig = figure(1)
ax = subplot(111)

data1 = loadtxt('problem1a1.dat',usecols=(0,2,3)) # t, x, v
data2 = loadtxt('problem1a2.dat',usecols=(0,2,3))
data3 = loadtxt('problem1a3.dat',usecols=(0,2,3))
x1= data1[:,1]
v1= data1[:,2]
x2= data2[:,1]
v2= data2[:,2]
x3= data3[:,1]
v3= data3[:,2]

xlim(-5,5)
ylim(-5,5)
trajectory1, = ax.plot([0],[0],'--', label="1,0")
trajectory2, = ax.plot([0],[0],'--', label="2,0")
trajectory3, = ax.plot([0],[0],'--', label="0,3")
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
show(block=False)

filenames = []

for i in range(len(data1)):
 trajectory1.set_xdata( x1[0:i] )
 trajectory1.set_ydata( v1[0:i] )
 trajectory2.set_xdata( x2[0:i] )
 trajectory2.set_ydata( v2[0:i] )
 trajectory3.set_xdata( x3[0:i] )
 trajectory3.set_ydata( v3[0:i] )
 time_text.set_text(time_template % (data1[:,0][i]))
 plt.legend()
 plt.pause(1e-30)
 draw()
 filename = f'{i}.png'
 filenames.append(filename)
 savefig(filename)

# build gif
with imageio.get_writer('mygif.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

# Remove files
for filename in set(filenames):
    os.remove(filename)


