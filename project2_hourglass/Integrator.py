from numpy import nan
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class VerletIntegrator(object):
  def __init__(self,dt,force):
    self.dt = dt
    self.force = force
    self.dims = force.dims
    #self.dims = 2
  def __call__(self,x,v,t):
    return self.forward(x,v,t)
  def forward(self,x,v,t):
    f = self.force
    a = f(x,v,t)
    dt = self.dt 
    dx = (v + .5*a*dt)*dt
    xn = x + dx
    vntemp = v  # estimate future v with present one
    dv = .5*(f(xn,vntemp,t+dt) + a)*dt  
    return (dx,dv)
  def backward(self,x,v,t):
    f = self.force
    self.dt = -self.dt
    (dx,dy) = self.forward(x,v,t)
    self.dt = -self.dt
    return (dx,dy)
  @property
  def L(self):
    return self.force.L
  def updateL(self,newL):
    return self.force.updateL(newL)
