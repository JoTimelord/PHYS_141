import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# References
# https://gist.github.com/neale/e32b1f16a43bfdc0608f45a504df5a84
# https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
# https://riptutorial.com/matplotlib/example/23558/basic-animation-with-funcanimation
# https://www.bragitoff.com/2020/10/3d-trajectory-animated-using-matplotlib-python/

# IMPORT DATA
data1 = np.loadtxt('Problem1.dat',usecols=(2,3,0)) # x, v, t
data2 = np.loadtxt('Problem11.dat',usecols=(2,3,0))
data3 = np.loadtxt('Problem12.dat',usecols=(2,3,0))

# ANIMATION FUNCTION
def func(num, dataSet, line):
    # NOTE: there is no .set_data() for 3 dim data...
    line.set_data(dataSet[0:2, :num])
    line.set_3d_properties(dataSet[2, :num])
    return line


# THE DATA POINTS
t = data1[:,2]
numDataPoints = len(t)

# GET SOME MATPLOTLIB OBJECTS
fig = plt.figure()
ax = Axes3D(fig)

# NOTE: Can't pass empty arrays into 3d version of plot()
line1 = plt.plot(data1[0], data1[1], data1[2], lw=2, c='g', label="1,0")[0] # For line plot
line2 = plt.plot(data2[0], data2[1], data2[2], lw=2, c='r', label="2,0")[0] # For line plot
line3 = plt.plot(data3[0], data3[1], data3[2], lw=2, c='b', label="0,3")[0] # For line plot

# AXES PROPERTIES]
ax.set_xlim3d([-5, 5])
ax.set_ylim3d([-5, 5])
ax.set_zlim3d([0, 10])
ax.set_xlabel('X(t)')
ax.set_ylabel('V(t)')
ax.set_zlabel('time')
ax.set_title('Trajectory of particles')

# Creating the Animation object
line_ani1 = animation.FuncAnimation(fig, func, frames=numDataPoints*10, fargs=(data1,line1), interval=200, blit=False)
line_ani2 = animation.FuncAnimation(fig, func, frames=numDataPoints*10, fargs=(data2,line2), interval=200, blit=False)
line_ani3 = animation.FuncAnimation(fig, func, frames=numDataPoints*10, fargs=(data3,line3), interval=200, blit=False)
#line_ani.save(r'AnimationNew.mp4')

plt.legend()
plt.show()
