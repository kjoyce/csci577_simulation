from Container import Container
from Force import SledForces
from Integrator import VerletIntegrator
from DistanceMatrix import DistanceMatrix
from numpy import linspace,sqrt,mod,inf,sin,pi
from pylab import rand,randint
class ParticleInitialize(object):
  def __init__(self,n_sled=13,load=0):
######   Params   ################
    xlim = (0,30)
    ylim = (-1,6)
    n_floor = 100
    start_sled = 0
    pull_force_lim = (-30,30)
    ave_vel_lim = (0,20)
    r = 2**(1./6.)
    dims = 2
    dt = .01
##################################
    
    self.case = "sled{}_load{}".format(n_sled,load)

    dy = r*sqrt(3)
    dx = r
    xinit = (2*start_sled+1)*.5*dx + (n_sled-1)*dx
    self.distance_matrix = DistanceMatrix()
    self.force = SledForces(dims,self.distance_matrix,xinit,n_sled,n_floor,float(load))
    self.integrate = VerletIntegrator(dt,self.force)
    self.c = Container(self.integrate,n_sled,n_floor)
    c = self.c

    for i in range(n_floor):
      c.addParticle(float(i)*dx,0.,0.,0.,1.)

    for i in range(n_sled):
      x = (2*start_sled+1)*.5*dx + i*dx
      y = dy*(i%2) + r
      c.addParticle(x,y,0,0.,1.)
      print "{}, {}".format(x,y)

    self.xlim            = xlim
    self.ylim            = ylim
    self.pull_force_lim  = pull_force_lim
    self.ave_vel_lim     = ave_vel_lim
  def __call__(self):
    return self.c,self.distance_matrix,self.force,self.integrate,self.xlim,self.ylim,self.pull_force_lim,self.ave_vel_lim,self.case


