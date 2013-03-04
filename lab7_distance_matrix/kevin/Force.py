from DistanceMatrix import PeroidicDistanceMatrix 
from numpy import array
class LeonardJonesForce(object):
  def __init__(self):
    self.sigma = 1
    self.eps = 1

  def __call__(self,container):
    Lx,Ly,Lz = container.Lx,container.Ly,container.Lz
    sigma = self.sigma
    eps = self.eps 
    distance_matrix = PeriodicDistanceMatrix(array([Lx,Ly,Lz]))
    dx = distance_matrix(array([container.xpos,container.ypos,container.zpos])
    r = abs(dx)
    return -14*eps/r*(2*(sigma/r)

