from pylab import pi,sqrt
from numpy import array
# Notice that parameters are global.
GM = 4*pi**2
m1 = 1
years = 100.
num_samples = 1000

# Initial Conditions for Planet
#Initial test, highly unstable orbit (falls straight into the star)
#Position
#Position
(x5,y5) = (3.6,0.)
#Velocity
(v5,w5) = (0, 6.4)
xinit = array([x5,y5,v5,w5]) 

# Saving filenames
figure_filename = 'medea_orbit5.png'
movie_filename = 'medea_orbit5.avi'

# Animation Params
frameskip = 1
time_length = 100.

