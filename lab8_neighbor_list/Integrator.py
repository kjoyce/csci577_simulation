from numpy import nan
from Force import LeonardJonesForce
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class VerletIntegrator(object):
  def __init__(self,dt,force):
    self.dt = dt
    self.force = force
    self.dims = force.dims
  def __call__(self,x,v,t):
    return self.forward(x,v,t)
  def forward(self,x,v,t):
    f = self.force
    a = f(x,v,t)
    dt = self.dt 
    dx = (v + .5*a*dt)*dt
    xn = x + dx
    vntemp = nan ### For a force that is velocity dependent, I don't know what to do, and this will break
    dv = .5*(f(xn,vntemp,t+dt,calc_auxilary=False) + a)*dt  
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

#### Unit Testing ####
if __name__ == '__main__':
  pass
