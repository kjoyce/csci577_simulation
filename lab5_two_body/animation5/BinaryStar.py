# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 18:17:33 2013

@author: Kevin Joyce, Nathan Sponberg
"""

from scipy.integrate import odeint # for integrate.odeint
from pylab import * # for plotting commands
from IPython.core.debugger import Tracer
debug_here = Tracer()
from matplotlib.pyplot import *
from matplotlib.mlab import find
from ODE_integrator import RungeKutta
from matplotlib import pyplot as plt
from matplotlib import animation
from binarystarparams import *
     
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
savefig(figure_filename)

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
                               frames=1000/frameskip, interval=time_length,
                               blit=True)

                  
plt.show()
anim.save(movie_filename,fps=15)#,extra_args=['-vcodec', 'libx264'])

