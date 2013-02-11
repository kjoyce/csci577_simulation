# -*- coding: utf-8 -*-
"""
Created on Sat Feb 09 16:03:27 2013
 
@author: Nathan Sponberg
"""
 
from scipy import *
from numpy import *
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib.mlab import find
debug_here = Tracer()
 
# Function performs Runga kutta method on another function passed to it
def runge_kutta(f,t,x,dt):
    k1 = f(t,x)*dt
    k2 = f(t+dt/2., x+k1/2.)*dt
    k3 = f(t+dt/2., x+k2/2.)*dt
    k4 = f(t+dt, x+k3)*dt
    return x + 1/6.*(k1+2*k2+2*k3+k4)     
 
def integrate(method, f, y0, t):
    y = ones([size(t), size(y0)])    
    y[0] = y[0]*y0
    for i in range(size(t)-1):
        y[i+1] = method(f,t[i],y[i],t[i+1]-t[i])
    return y
 
def FallingBody(t,x):
    return array([x[1], -9.8])
 
def fallWithDrag(t,x):
    ratio = x[1]/terminalv
    return array([x[1], -9.8*(1-ratio)])
 
def fallWithDragSquared(t,x):
    ratio = x[1]/terminalv**2
    return array([x[1], -9.8*(1-ratio)])
 
# Time in seconds: data[0.:]
# y-axis position in meters: data[1,:]
 
data = array([[.2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,.4280,
              .4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,.6498,.6744,
              .6990,.7236,.7482,.7728,.7974,.8220,.8466],
              [.4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,.3609,
               .3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,.2497,.2337,
               .2175,.2008,.1846,.1696,.1566,.1393,.1263]])
 
################Simulation##############
 
terminalv = -.4
x0 = data[1][0]
xprime0 = (data[1][1] - data[1][0])/(data[0][1]-data[0][0])
(a,b) = (data[0][0],data[0][-1])
startState = [array([x0,xprime0])]        
iterations = 100
dt = (b-a)/iterations
t = linspace(a,b,iterations)
 
fallingNumericSolution = integrate(runge_kutta, fallWithDrag, startState, t)
fallingNumericSolution2 = integrate(runge_kutta, fallWithDragSquared, startState, t)
idx = find(fallingNumericSolution2[:,0] > 0)
fallingNumericSolution2 = fallingNumericSolution2[idx]
 
figure(1)
plot(t,fallingNumericSolution[:,0])
plot(t[idx],fallingNumericSolution2[:,0])
plot(data[0],data[1])
title("Falling Bodies With Drag")
xlabel("Time (sec)")
ylabel("Displacement (meters)")
legend(("Linear", "Quadratic", "Emperical Data"))
show()
