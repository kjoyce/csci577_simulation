from pylab import pi,sqrt
from numpy import array
# Notice that parameters are global.
GM = 4*pi**2
m1 = 1
years = 100.
num_samples = 1000

# Even more stable eliptical orbit farther away from binary system,
#starts farther out from the stars (4 AUs on the y axis)
#Position
(x4,y4) = (1.1,4.)
#Velocity
(v4,w4) = (-4.,2)
xinit = array([x4,y4,v4,w4]) 

# Saving filenames
figure_filename = 'orbit4.png'
movie_filename = 'orbit4.avi'

# Animation Params
frameskip = 1
time_length = 100.

