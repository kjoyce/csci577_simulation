from pylab import pi,sqrt
from numpy import array
# Notice that parameters are global.
GM = 4*pi**2
m1 = 1
years = 100.
num_samples = 1000

#Stable eliptical orbit passes very close to binary system,
#starts farther out from the stars (2.1 AUs on the x axis)

#Position
(x3,y3) = (2.1,1.)
#Velocity
(v3,w3) = (-9.,3)
xinit = array([x3,y3,v3,w3]) 

# Saving filenames
figure_filename = 'orbit3.png'
movie_filename = 'orbit3.avi'

# Animation Params
frameskip = 1
time_length = 100.

