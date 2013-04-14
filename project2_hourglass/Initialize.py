from Container import Container
from Force import HourGlassForce
from Integrator import VerletIntegrator
from DistanceMatrix import DistanceMatrix
from numpy import linspace,sqrt,mod,inf,sin,pi
from pylab import rand,randint
class ParticleInitialize(object):
  def __init__(self):
######   Params   ################
    xlim = (-5,5)
    ylim = (-1,8)
    n_floor = 10
    r = 2**(1./6.)
    dims = 2
    dt = .01
##################################
    
    self.case = "Test case"

    self.distance_matrix = DistanceMatrix()
    self.force = HourGlassForce(self.distance_matrix)
    self.integrate = VerletIntegrator(dt,self.force)
    self.c = Container(self.integrate,n_floor)

    for x in linspace(xlim[0],xlim[1],n_floor):
      print x
      self.c.addParticle(x,0.,0.,0.,1.)

    self.c.addParticle(0.,1.5,0.,0.,1.) 
    self.xlim            = xlim
    self.ylim            = ylim

  def __call__(self):
    return self.c,self.distance_matrix,self.force,self.integrate,self.xlim,self.ylim,self.case


