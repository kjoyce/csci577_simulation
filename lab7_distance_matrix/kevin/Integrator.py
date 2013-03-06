from numpy import nan
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class VerletIntegrator(object):
  def __init__(self,dt):
    self.dt = dt

  def __call__(self,f,x,v,t):
    return self.forward(f,x,v,t)

  def forward(self,f,x,v,t):
    dt = self.dt
    dx = v*dt + .5 * f(x,v,t)*dt**2
    xn = x + dx
    vntemp = nan ### For a force that is velocity dependent, I don't know what to do, and this will break
    dv = .5*(f(xn,vntemp,t+dt) + f(x,v,t))*dt  
    return (dx,dv)

  def backwards(self,f,x,v,t):
    pass

#### Unit Testing ####
if __name__ == '__main__':
  pass
