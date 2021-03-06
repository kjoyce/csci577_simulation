from scipy import *
from numpy import *
from IPython.core.debugger import Tracer
debug_here = Tracer()
def numerical_integrate(method,f,y0,t):
  # maybe make this different in the case that y0 is 1-dimensional
  y = ones([size(t),size(y0)])
  y[0] = y[0]*y0 
  for i in range(size(t)-1):
      y[i+1] = method(t,y[i],f,(t[i+1] - t[i]))
  return y

def euler(t,y,f,dt):
  return  y + f(t,y)*dt

# 1. Falling Object
(a,b) = (0,1)
def falling_object_system(deriv,m):
  return array([ deriv[1], -9.8/m ])

def numeric_falling_object(t,m):
  return numerical_integrate(euler,lambda t,y: falling_object_system(y,m),

def analytic_falling_object(t,m):
  return #

# 2. Harmonic Oscillator
def oscilator_system(deriv,k):
  return #
  
