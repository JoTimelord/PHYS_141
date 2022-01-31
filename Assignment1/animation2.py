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

data1 = loadtxt('planetorbitdat/Earth.dat',usecols=(0,1,3)) #t, x, y
t= data1[:,0]
x= data1[:,1]
y= data1[:,2]

xlim(-2,2)
ylim(-2,2)
trajectory, = ax.plot([0],[0],'--', label="earth's orbit")
planet,      = ax.plot([0],[0],'or', ms=10)
sun,         = ax.plot([0],[0], 'oy', ms=50)
time_template = 'time = %.1fyear'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
show(block=False)

filenames = []

for i in range(len(data1)):
 planet.set_xdata(x[i])
 planet.set_ydata(y[i])
 trajectory.set_xdata( x[0:i] )
 trajectory.set_ydata( y[0:i] )
 time_text.set_text(time_template % (data1[:,0][i]))
 plt.pause(1e-30)
 draw()
 filename = f'{i}.png'
 filenames.append(filename)
 savefig(filename)

# build gif
with imageio.get_writer('earthanimation.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

# Remove files
for filename in set(filenames):
    os.remove(filename)


