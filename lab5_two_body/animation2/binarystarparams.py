from pylab import pi,sqrt
from numpy import array
# Notice that parameters are global.
GM = 4*pi**2
m1 = 1
years = 100.
num_samples = 1000

#Changing the velocity so that it is in the right direction,
#somewhat more stable orbit passing very close to stars
#Position
(x2,y2) = (1.1,1.)
#Velocity
(v2,w2) = (-8.,0.2)
xinit = array([x2,y2,v2,w2]) 

# Saving filenames
figure_filename = 'orbit2.png'
movie_filename = 'orbit2.avi'

# Animation Params
frameskip = 1
time_length = 100.

