from numpy import nan
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class VerletIntegrator(object):
  def __init__(self,dt,force):
    self.dt = dt
    self.force = force

  def __call__(self,x,v,t):
    return self.forward(f,x,v,t)

  def forward(self,x,v,t):
    f = self.force
    a = f(x,v,t)
    dt = self.dt 
    dx = (v + .5*a*dt)*dt
    xn = x + dx
    vntemp = nan ### For a force that is velocity dependent, I don't know what to do, and this will break
    dv = .5*(f(xn,vntemp,t+dt,calc_auxilary=False) + a)*dt  
    return (dx,dv)

  # I am not sure about this..
  def backward(self,x,v,t):
    f = self.force
    self.dt = -self.dt
    (dx,dy) = self.forward(f,x,v,t)
    self.dt = -self.dt
    return (dx,dy)


#### Unit Testing ####
if __name__ == '__main__':
  pass
