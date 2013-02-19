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
(x1,y1) = (1.1,1.)
#Velocity
(v1,w1) = (sqrt(GM/y1),sqrt(GM/x1))
xinit = array([x1,y1,v1,w1]) 

# Saving filenames
figure_filename = 'orbit1.png'
movie_filename = 'orbit1.avi'

# Animation Params
frameskip = 1
time_length = 100.

