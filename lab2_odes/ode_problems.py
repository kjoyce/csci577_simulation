# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 15:49:39 2013

@author: Nathan Sponberg, Kevin Joyce, Patrick Funk
"""

from scipy import *
from numpy import *
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
def FallingBody(x,dt):
    return array([x[1], -9.8])

# Function simulates simple harmonic oscillation, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
def Oscillator(x,dt):
    k=1
    return array([x[1], -k*x[0]])

# Analytic solution for a falling body
def FallingAnalytic(t,y0):
  yprime0 = 0 # hack for now
  return array([-1./2. * 9.8 * t**2 + y0,-9.8*t+yprime0]).T
  #return -1./2. * 9.8 * t**2 + y0

# Analyitic solution for Oscillation

def OscilAnalytic(k, m, x0, dt, numPoints):
    #result = []    
    #for i in range(numPoints):
    #    result.append([x0*math.cos(math.sqrt(k/m)*(i*dt))])
    #return result
    ## This version takes advantage of the numpy array
    t = linspace(0,dt*numPoints,numPoints)
    return array([x0*cos(sqrt(k/m)*t), -x0*sqrt(k/m)*sin(sqrt(k/m)*t)]).T


def OneDMotionSim(f, startState, resolution, dt):
    states = startState
    for i in range(resolution):
        states.append(euler(states[-1], f, dt))
    return states

################# Begin Simulations ######################
# 1 Falling object 
#initial state for falling object
fall_resolution = 1000
(a,b) = (0.,2.)
dt = (b-a)/fall_resolution
stateFalling = [array([10.0,0.0])]
fall_t = linspace(a,b,fall_resolution)
fall_analytic_solution = FallingAnalytic(fall_t,stateFalling[0][0])
numeric_states = OneDMotionSim(FallingBody,stateFalling,fall_resolution-1,dt)
fall_numeric_solution = array(numeric_states)
# obtain only positive values

# Total Energy for simulation and analytic solution
E_fall_simulation = 9.8*fall_numeric_solution[:,0] + .5*fall_numeric_solution[:,1]**2
E_fall_analytic = 9.8*fall_analytic_solution[:,0] + .5*fall_analytic_solution[:,1]**2
# Percent difference
E_fall_change = abs((E_fall_simulation[0]-E_fall_simulation[-1])/E_fall_simulation[0])

# plot velocity and time
idx = find( fall_analytic_solution[:,0] >= 0 )
figure()
plot(fall_t[idx],fall_analytic_solution[idx,0])
plot(fall_t[idx],fall_numeric_solution[idx,0])
title("Falling Object, n={:d}".format(fall_resolution))
xlabel("Time (sec)")
ylabel("Height Above Ground (meters)")
legend(("Analytic","Numeric"))

figure()
plot(fall_t[idx],fall_analytic_solution[idx,1])
plot(fall_t[idx],fall_numeric_solution[idx,1])
title("Velocity of Object, n={:d}".format(fall_resolution))
xlabel("Time (sec)")
ylabel("Velocity (meters/sec)")
legend(("Analytic","Numeric"))

# plot of Energy vs. time
figure()
plot(fall_t,E_fall_simulation)
plot(fall_t,E_fall_analytic)
title("Energy vs. time of Falling Obect, n={:d}".format(fall_resolution))
xlabel("Time (sec)")
ylabel("Energy (J)")
legend(("Analytic","Numeric"))
comment = r"$ \left|\frac{{E(0) - E(t_{{end}})}}{{E(0)}}\right| \approx {0:.4f} $".format(E_fall_change)
text(0,98.175,comment,fontsize=20)
show()

################################
# 2. Simple Harmonic Motion
#initial state for Oscillation
k = 1.
m = 1.
x0 = 1.
xprime0 = 0.
(a,b) = (0.,40.)
startOscil = [array([x0,xprime0])]

#runs simulation for oscillator, time interval is 40000*0.001 = 40 secs
osc_resolution = 40000
dt = (b-a)/osc_resolution
oscilStates = OneDMotionSim(Oscillator, startOscil, osc_resolution-1, dt)
num_oscil_solution = array(oscilStates)
osc_t = linspace(a,b,osc_resolution)
#note: dt and numPoints should match dt and osc_resolution
#from the simulation
osc_analytic = OscilAnalytic(k,m,x0,dt,osc_resolution);

# Total Energy for spring

E_spring_total = 1./2.*k*x0

# Total Energy for simulation

E_spring_simulation = 1./2.*k*x0*array(num_oscil_solution[:,0])**2 + 1./2.*m*array(num_oscil_solution[:,1])**2
E_spring_analytic = 1./2.*k*x0*array(osc_analytic[:,0])**2 + 1./2.*array(osc_analytic[:,1])**2


# percent difference
E_spring_change = abs(E_spring_total-E_spring_simulation[-1])/(E_spring_total)


figure()
plot(osc_t,num_oscil_solution[:,0])
plot(osc_t,osc_analytic[:,0])
title("Simple Harmonic Motion, n={:d}".format(osc_resolution))
xlabel("Time (sec)")
ylabel("Displacement (meters)")
legend(("Analytic","Numeric"))

figure()
plot(osc_t,num_oscil_solution[:,1])
plot(osc_t,osc_analytic[:,1])
title("Velocity of Simple Harmonic Motion, n={:d}".format(osc_resolution))
xlabel("Time (sec)")
ylabel("Velocity (meters/sec)")
legend(("Analytic","Numeric"))

# plot of Energy
figure()
plot(osc_t,E_spring_simulation)
plot(osc_t,E_spring_analytic)
title("Energy vs. Time for Spring, n={:d}".format(osc_resolution))
legend(("Analytic","Numeric"))
comment = r"$ \left|\frac{{E(0) - E(t_{{end}})}}{{E(0)}}\right| \approx {0:.4f} $".format(E_spring_change)
text(7,.52,comment,fontsize=20)

show()
