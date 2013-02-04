# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 15:49:39 2013

@author: Nathan Sponberg
"""

from scipy import *
from numpy import *
import pylab as pl
from IPython.core.debugger import Tracer
debug_here = Tracer()

# Function performs eulers method on another passed to it
def euler(x, f, dt):
    return x+f(x,dt)*dt

# Function simulates falling object, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
def FallingBody(x,dt):
    return array([(x[1]-(9.8*dt)), -9.8])

# Function simulates simple harmonic oscillation, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
def Oscillator(x,dt):
    k=1
    return array([x[1], -k*x[0]])

# Analytic solution for a falling body
# Broken! object falls too fast, needs to be fixed
#def FallingAnalytic(a, y0, dt, numPoints):
#    result = [y0]
#    for i in range(numPoints):
#        result.append(result[-1]-0.5*a*(dt*(i+1))**2)
#    return result

# Analyitic solution for Oscillation

def OscilAnalytic(k, m, x0, dt, numPoints):
    result = []    
    for i in range(numPoints):
        result.append([x0*math.cos(math.sqrt(k/m)*(i*dt))])
    return result


#initial state for falling object
stateFalling = [array([10.0,0.0])]

# simulates falling object until it hits the ground
#while stateFalling[-1][0] > 0:
#    stateFalling.append(euler(stateFalling[-1], FallingBody, 0.1))
    


def OneDMotionSim(f, startState, resolution, dt):
    states = startState
    for i in range(resolution):
        states.append(euler(states[-1], f, dt))
    return states

#initial state for Oscillation
startOscil = [array([1.0,0.0])]

#runs simulation for oscillator, time interval is 40000*0.001 = 40 secs
oscilStates = OneDMotionSim(Oscillator, startOscil, 40000, 0.001)

#plots the oscillation with position and velocity
pl.plot(oscilStates)

#plots the analytic solution for the oscillation
#note: dt and numPoints should match dt and resolution
#from the simulation
pl.plot(OscilAnalytic(1,1,1,0.001,40000))

pl.show()
