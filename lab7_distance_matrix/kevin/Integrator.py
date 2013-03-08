from numpy import nan
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class VerletIntegrator(object):
  def __init__(self,dt):
    self.dt = dt

  def __call__(self,f,x,v,t):
    return self.forward(f,x,v,t)

  def forward(self,f,x,v,t):
    a = f(x,v,t)
    dt = self.dt 
    dx = (v + .5*a*dt)*dt
    xn = x + dx
    vntemp = nan ### For a force that is velocity dependent, I don't know what to do, and this will break
    dv = .5*(f(xn,vntemp,t+dt) + a)*dt  
    return (dx,dv)

  # I am not sure about this..
  def backward(self,f,x,v,t):
    #self.dt = -self.dt
    #(dx,dy) = self.forward(f,x,v,t)
    (dx,dy) = self.forward(f,x,v,t)
    (dx,dy) =  (-dx,-dy)
    #fback = lambda x,v,t: -f(x,v,t) 
    #(dx,dy) = self.forward(fback,x,v,t)
    return (dx,dy)


#### Unit Testing ####
if __name__ == '__main__':
  pass
