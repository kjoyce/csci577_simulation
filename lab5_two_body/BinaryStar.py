# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 18:17:33 2013

@author: Kevin Joyce, Nathan Sponberg
"""

from scipy.integrate import odeint # for integrate.odeint
from pylab import * # for plotting commands
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
from ODE_integrator import RungeKutta
from matplotlib import pyplot as plt
from matplotlib import animation
debug_here = Tracer()
 
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
xinit = InitStates[2]            

#Two Body function with the second body fixed (second star)
def PlanetOrbit(state,t): #OUCH! The signature is reversed for odeint!
    x1 = array([state[0],state[1]]) # body 1 position
    v1 = array([state[2],state[3]]) # body 1 velocity
    star2 = array([2.0,0.]) #star number 2

    r1 = norm(x1)
    r21 = norm(star2 - x1)

    a1 = -GM/r1**3*x1 + m1*GM/r21**3*(star2-x1)
    derivatives = array([v1,a1]).T.flatten(1)
    return derivatives
 
times = linspace(0.0,years,num_samples)
scipy_result = odeint(PlanetOrbit,xinit,times)
scipy_result = scipy_result.T

#######################Plot and Animation##########################

#Static plot
figure(figsize = (6,6))
plot(0,0,'ro')
plot(2,0,'yo')
plot(scipy_result[0,:],scipy_result[1,:])
xlim((-9,11))
ylim((-10,10))
grid()
title('Planetary Orbit')
xlabel('Horizontal Distance (AU)')
ylabel('Vertical Distance (AU)')
comment = r"""Initial Conditions
$x_0 = {0}$
$y_0 = {1}$
$x_{{v0}} = {2}$
$y_{{v0}} = {3}$""".format(xinit[0],xinit[1],xinit[2],xinit[3])
text(-8,6,comment,fontsize=10)

#Animation
fig = plt.figure(figsize = (6,6))
ax = plt.axes(xlim=(-9, 11), ylim=(-10, 10))
line, = ax.plot([], [])
point, = ax.plot([], [], 'g^')
ax.plot(0,0,'ro')
ax.plot(2,0,'yo')
title('Planetary Orbit Animation')
xlabel('Horizontal Distance (AU)')
ylabel('Vertical Distance (AU)')

#Number of frames to skip (n-1)
frameskip = 1

#Set to true to trace planetary orbit in animation
trace = False

#Initial animation state
def init():
    line.set_data([], [])
    point.set_data([], [])
    return line, point,

#Animation loop
def animate(i):
    saveImg = True
    x = scipy_result[0,i]
    y = scipy_result[1,i]
    line.set_data(x, y)
    return line,
    
def animateTrace(i):
    if trace:
        x = scipy_result[0,0:frameskip*i]
        y = scipy_result[1,0:frameskip*i]
        line.set_data(x, y)
        
    px = scipy_result[0,frameskip*i]
    py = scipy_result[1,frameskip*i]    
    point.set_data(px,py)
    return line, point,
    
#Start animation
anim = animation.FuncAnimation(fig, animateTrace, init_func=init,
                               frames=1000/frameskip, interval=40,
                               blit=True)
  
                  
plt.show()

