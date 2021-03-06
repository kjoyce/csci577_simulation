# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 15:49:39 2013

@author: Nathan Sponberg
"""

from scipy import *
from numpy import *
import pylab as pl
from matplotlib import *
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
def FallingAnalytic(t,y0):
    return -1/2*-9.8*t+y0

# Analyitic solution for Oscillation

def OscilAnalytic(k, m, x0, dt, numPoints):
    result = []    
    for i in range(numPoints):
        result.append([x0*math.cos(math.sqrt(k/m)*(i*dt))])
    return result


#initial state for falling object
stateFalling = [array([10.0,0.0])]   


def OneDMotionSim(f, startState, resolution, dt):
    states = startState
    for i in range(resolution):
        states.append(euler(states[-1], f, dt))
    return states

#########  Simulations  ###########

#initial state for falling object
stateFalling = [array([10.0,0.0])]
   

startFall = OneDMotionSim(FallingBody, stateFalling, 40000, 0.001)

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

#####   Energy    ##########

# Total Energy for Falling object

E_fall_total = 9.8*10

# Total Energy for simulation

E_fall_simulation = 9.8*array(startFall[:,0]) + .5*array(startFall[:,1])**2

# plot of result
figure(1)
pl.plot(E_fall_simulation)
axhline(y=98)

# Percent difference
E_fall_change = abs(E_fall_total-E_fall_simulation(end))/(E_fall_total)

# Total Energy for spring

E_spring_total = 1/2

# Total Energy for simulation

E_spring_simulation = 1/2*array(oscilStates[:,0])**2 + 1/2*array(oscilStates[:,1])

# plot of result
figure(2)
pl.plot(E_spring_simulation)
axhline(y=1/2)

# percent difference

E_spring_change = abs(E_spring_total-E_spring_simulation(end))/(E_spring_total)
