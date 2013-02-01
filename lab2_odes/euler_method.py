from scipy import *
from numpy import *
from IPython.core.debugger import Tracer
debug_here = Tracer()
def numerical_integrate(f,y0,a,b,dt):    
  t = arange(a,b,dt)    
  y = ones([size(t),size(y0)])
  y[0] = y[0]*y0 
  for i in range(size(t)-1):
      y[i+1] = euler(y[i],f,dt)
  return (t,y)

def euler(y,f,dt):
  return  #y + f(t,y)*dt

def newtons_second_law(t,y,F,m):
  return array([ y[1], F/m ])

def analytic_solution(t):
  return
