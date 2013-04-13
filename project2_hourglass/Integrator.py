from numpy import nan
from Force import LeonardJonesForce
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
### Test integrating a spring between two particles ##
  from numpy import array,sqrt,zeros,arange
  from DistanceMatrix import DistanceMatrix
  from pylab import figure,plot,show
  distance_matrix = DistanceMatrix()
  x = array([[0,0],[sqrt(2)**-1,sqrt(2)**-1]])
  v = array([[0,0],[0,0]])
  dx = distance_matrix(x)
  print x
  print dx.T
  def a(x,v,t,calc_auxilary=False):
    return 3*x
  t = arange(0,4,.1)
  xx = zeros(t.shape)
  yy = zeros(t.shape)
  xxx = zeros(t.shape)
  yyy = zeros(t.shape)
  integrate = VerletIntegrator(.1,a)
  for i in range(len(t)):
    xx[i],yy[i] = x[0,0],x[0,1]
    xxx[i],yyy[i] = x[1,0],x[1,1]
    (dx,dv) = integrate(x,v,i)
    x += dx
    v += dv
  figure()
  plot(t,xx)
  plot(t,xxx)
  plot(t,yy)
  plot(t,yyy)
  show()
