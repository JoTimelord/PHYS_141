#!/bin/bash

set PATH = /usr/local/physics/bin/glnemo2.app/Contents/MacOS:$PATH

csh

magalie gal.ic ndsik=10000 nbulge=0 nhalo=1000
whoami
# The first ten thousand particles (indices 0:9999) are the disk:
glnemo2 gal.ic 0:9999 xrot=90
# The next ten thousand particles (indices 10000:19999) are the halo:
glnemo2 gal.ic 10000:19999 xrot=90
# run the simulation
gyrfalcON gal.ic gal.out tstop=50 step=0.1 kmax=6 eps=0.03

glnemo2 gal.out 0:9999,10000:19999
mkdir frames
cd frames

glnemo2 ../gal.out 0:9999 hsize=400 wsize=400 play=t screenshot=frame
ffmpeg -y -i frame.%05d.jpg -c:v libx264 -vf fps=20 -pix_fmt yuv420p
video.mp4
