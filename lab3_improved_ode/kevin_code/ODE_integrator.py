from numpy import arange,array
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class ODE_integrator(object):
  """ A class for integrating ODE's of the form dx/dt = f(x,t) over a finite interval (a,b).  By default, this performs fourth order Runge-Kutta. """
  def __init__(self,f,y0,a,b,dt=.1):
    self.f = f
    self.a = a
    self.b = b
    self.y = [y0]
    self.t = [a]
    self.dt = dt

  def integrate(self):
    while( self.t[-1] < self.b ): 
      (t,y) = self.update_function()
      self.t.append(t)
      self.y.append(y)
    self.y = array(self.y).T
    return (array(self.t),self.y)

  def update_function(self):
    (f,y,t,dt) = (self.f,self.y[-1],self.t[-1],self.dt)
    k1 = f(t,y)*dt
    k2 = f(t + dt/2,y + k1/2)*dt
    k3 = f(t + dt/2,y + k2/2)*dt
    k4 = f(t + dt,y + k3)*dt
    return (t+dt,y + 1./6.*(k1 + 2*k2 + 2*k3 + k4))

class RungeKutta(ODE_integrator):
  pass # This is the default behavior, so do nothing
  
class Euler(ODE_integrator):
  def update_function(self):
    (f,y,t,dt) = (self.f,self.y[-1],self.t[-1],self.dt)
    return (t+dt,y+f(t,y)*dt)

class EulerRichardson(ODE_integrator):
  def update_function(self):
    (f,y,t,dt) = (self.f,self.y[-1],self.t[-1],self.dt)
    return (t+dt,y + f(t+dt/2,y+ f(t,y)*dt/2)*dt)

class PredictorCorrector(ODE_integrator):
  def __init__(self,f,y0,a,b,dt='.1'):
    super(PredictorCorrector, self).__init__(f,y0,a,b,dt)
# spin up with RungeKutta
    (t,y) = super(PredictorCorrector, self).update_function()
    self.t.append(t)
    self.y.append(y)
  def update_function(self):
    (f,y,t,dt) = (self.f,self.y[-1],self.t[-1],self.dt)
    yp = self.y[-2] + 2*f(t,y)*dt
    return (t+dt,y + .5*(f(t,yp) + f(t,y))*dt)
