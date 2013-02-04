# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 15:49:39 2013

@author: Nathan Sponberg
"""

from scipy import *
from numpy import *
import pylab as pl
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
debug_here = Tracer()

# Function performs eulers method on another passed to it
def euler(x, f, dt):
    return x+f(x,dt)*dt

# Function simulates falling object, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
# Function simulates simple harmonic oscillation, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
def Oscillator(x,dt):
    k=1
    return array([x[1], -k*x[0]])

def FallingBody(x,dt):
    return array([x[1], -9.8])

# Analytic solution for a falling body
def FallingAnalytic(t,y0):
  return -1./2. * 9.8 * t**2 + y0

# Analyitic solution for Oscillation

def OscilAnalytic(k, m, x0, dt, numPoints):
    result = []    
    for i in range(numPoints):
        result.append([x0*math.cos(math.sqrt(k/m)*(i*dt))])
    return result


def OneDMotionSim(f, startState, resolution, dt):
    states = startState
    for i in range(resolution):
        states.append(euler(states[-1], f, dt))
    return states

################# Begin Simulations ######################
# 1 Falling object 
#initial state for falling object
resolution = 1000
stateFalling = [array([10.0,0.0])]
(a,b) = (0.,2.)
t = linspace(a,b,resolution)
x_exact = FallingAnalytic(t,stateFalling[0][0])
numeric_states = OneDMotionSim(FallingBody,stateFalling,resolution-1,(b-a)/resolution)
numeric_solution = array(numeric_states)

figure(1)
clf()
# obtain only positive values
idx = find( x_exact >= 0 )
plot(t[idx],x_exact[idx])
#debug_here()
plot(t[idx],numeric_solution.T[0][idx])
legend(("Analytic Solution","Numeric Solution"))

#initial state for Oscillation
startOscil = [array([1.0,0.0])]

#runs simulation for oscillator, time interval is 40000*0.001 = 40 secs
oscilStates = OneDMotionSim(Oscillator, startOscil, 40000, 0.001)

figure(2)
#plots the oscillation with position and velocity
pl.plot(oscilStates)

#plots the analytic solution for the oscillation
#note: dt and numPoints should match dt and resolution
#from the simulation
pl.plot(OscilAnalytic(1,1,1,0.001,40000))

pl.show()
