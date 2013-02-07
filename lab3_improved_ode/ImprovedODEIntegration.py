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
def euler(f,t,x,dt):
    return x+f(t,x)*dt

# Function performs euler-richardson method on another passed to it
def euler_richardson(f,t,x,dt):
    return x + f(t+dt/2,x+f(t,x)*dt/2)*dt
    
# Function performs Runga kutta method on another function passed to it
def runge_kutta(f,t,x,dt):
    k1 = f(t,x)*dt
    k2 = f(t+dt/2., x+k1/2.)*dt
    k3 = f(t+dt/2., x+k2/2.)*dt
    k4 = f(t+dt, x+k3)*dt
    return x + 1/6.*(k1+2*k2+2*k3+k4)        
    
def predictor_corrector(f,t,x,previousState,dt):
    xPredict = previousState + 2.*f(t,x)*dt
    return x + 1./2.*(f(t,xPredict) + f(t,x))*dt
    

# Function simulates simple harmonic oscillation, x is a vector contain
# position and velocity data, outputs velocity and acceleration
# dt is the time step
def Oscillator(t,x):
    k=1
    return array([x[1], -k*x[0]])

def OscilAnalytic(k, m, x0, dt, numPoints):
    #result = []
    #for i in range(numPoints):
    # result.append([x0*math.cos(math.sqrt(k/m)*(i*dt))])
    #return result
    ## This version takes advantage of the numpy array
    t = linspace(0,dt*numPoints,numPoints)
    return array([x0*cos(sqrt(k/m)*t), -x0*sqrt(k/m)*sin(sqrt(k/m)*t)]).T

def oscil_energy_calc(k,m,x0,states):
    return 1./2.*k*x0*array(states[:,0])**2 + 1./2.*m*array(states[:,1])**2
    
def energy_percent_change(simEnergy,totalEnergy):
    return abs(totalEnergy-simEnergy[-1])/(totalEnergy)
    
def integrate(method, f, y0, t):
    y = ones([size(t), size(y0)])    
    if method == predictor_corrector:
        y[0] = y[0]*y0
        y[1] = runge_kutta(f,t,y[0],t[1] - t[0])
        for i in range(1,size(t)-2):
            y[i+1] = method(f,t[i],y[i],y[i-1],t[i+1]-t[i])
        return y
    else :
        y[0] = y[0]*y0
        for i in range(size(t)-1):
            y[i+1] = method(f,t[i],y[i],t[i+1]-t[i])
        return y

################################
# 2. Simple Harmonic Motion
#initial state for Oscillation
k = 1.
m = 1.
x0 = 1.
xprime0 = 0.
(a,b) = (0.,40.)
startOscil = [array([x0,xprime0])]

#runs simulation for oscillator, dt is resolution divided by total 
#time interval
osc_resolution = 270
dt = (b-a)/osc_resolution
t = linspace(a, b, osc_resolution)
num_oscil_solution1 = integrate(euler, Oscillator, startOscil, t)
num_oscil_solution2 = integrate(euler_richardson, Oscillator, startOscil, t)
num_oscil_solution3 = integrate(runge_kutta, Oscillator, startOscil, t)
num_oscil_solution4 = integrate(predictor_corrector, Oscillator, startOscil, t)

#note: dt and numPoints should match dt and osc_resolution
#from the simulation
osc_analytic = OscilAnalytic(k,m,x0,dt,osc_resolution);

# Total Energy for spring

E_spring_total = 1./2.*k*x0

# Total Energy for simulation

eulerEnergy = oscil_energy_calc(k,m,x0,num_oscil_solution1)
eulerRichEnergy = oscil_energy_calc(k,m,x0,num_oscil_solution2)
rungeKuttaEnergy = oscil_energy_calc(k,m,x0,num_oscil_solution3)
predictorEnergy = oscil_energy_calc(k,m,x0,num_oscil_solution4)
analyticEnergy = oscil_energy_calc(k,m,x0,osc_analytic)

eulerPercentChange = energy_percent_change(eulerEnergy,E_spring_total)
eulerRichPercentChange = energy_percent_change(eulerRichEnergy,E_spring_total)
rungePercentChange = energy_percent_change(rungeKuttaEnergy,E_spring_total)
predictorPercentChange = energy_percent_change(predictorEnergy,E_spring_total)

#E_spring_simulation = 1./2.*k*x0*array(num_oscil_solution[:,0])**2 + 1./2.*m*array(num_oscil_solution[:,1])**2
#E_spring_analytic = 1./2.*k*x0*array(osc_analytic[:,0])**2 + 1./2.*array(osc_analytic[:,1])**2


# percent difference
#E_spring_change = abs(E_spring_total-E_spring_simulation[-1])/(E_spring_total)


figure(1)
#plot(t,num_oscil_solution1[:,0])
plot(t,num_oscil_solution2[:,0])
#plot(t,num_oscil_solution3[:,0])
#plot(t,num_oscil_solution4[:,0])
plot(t,osc_analytic[:,0])
title("Simple Harmonic Motion, n={:d}".format(osc_resolution))
xlabel("Time (sec)")
ylabel("Displacement (meters)")
legend(("Euler", "Euler Richardson", "Runge Kutta", "Analytic"))
#show()
'''
figure()
plot(osc_t,num_oscil_solution[:,1])
plot(osc_t,osc_analytic[:,1])
title("Velocity of Simple Harmonic Motion, n={:d}".format(osc_resolution))
xlabel("Time (sec)")
ylabel("Velocity (meters/sec)")
legend(("Analytic","Numeric"))
'''

#plot of Energy
figure()
plot(t,analyticEnergy)
plot(t,eulerRichEnergy)
title("Energy vs. Time for Spring, n={:d}".format(osc_resolution))
legend(("Analytic","Numeric"))
comment = r"$ \left|\frac{{E(0) - E(t_{{end}})}}{{E(0)}}\right| \approx {0:.4f} $".format(eulerRichPercentChange)
text(7,.52,comment,fontsize=20)

show()



