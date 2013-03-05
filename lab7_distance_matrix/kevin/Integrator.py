class VerletIntegrator(object):
  def __init__(self,dt):
    self.dt = dt

  def __call__(self,f,x,v,t):
    dt = self.dt
    xn = x + v*dt + .5 * f(x,v,t)*dt**2 
    vntemp = nan ### For a force that is velocity dependent, I don't know what to do, and this will break
    vn = v + .5(f(xn,vntemp,t+dt) + f(x,v,t))*dt  

#### Unit Testing ####
if __name__ == '__main__':
  pass
