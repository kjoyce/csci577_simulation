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

#Changing the velocity so that it is in the right direction,
#somewhat more stable orbit passing very close to stars

#Position
(x2,y2) = (1.1,1.)
#Velocity
(v2,w2) = (-8.,0.2)

#Stable eliptical orbit passes very close to binary system,
#starts farther out from the stars (2.1 AUs on the x axis)

#Position
(x3,y3) = (2.1,1.)
#Velocity
(v3,w3) = (-9.,3)

# Even more stable eliptical orbit farther away from binary system,
#starts farther out from the stars (4 AUs on the y axis)

#Position
(x4,y4) = (1.1,4.)
#Velocity
(v4,w4) = (-4.,2)

#Very stable orbit, almost circular stays about 1.6 to 2 AUs from stars
#starts at 3.6 AUs along the x axis (this is 1.6 AUs from the yellow star on the
#x axis, this creates a stable orbit)

#Position
(x5,y5) = (3.6,0.)
#Velocity
(v5,w5) = (0, 6.4)

#Loads the initial states into a list
InitStates = [array([x1,y1,v1,w1]), array([x2,y2,v2,w2]),
              array([x3,y3,v3,w3]), array([x4,y4,v4,w4]),
              array([x5,y5,v5,w5])]
              
#Chose the initial state you want by changing the index            
xinit = InitStates[1]            
save_filename = 'orbit1.avi'

frameskip = 1

time_length = 100.

