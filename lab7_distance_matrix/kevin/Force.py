""" This module provides functions for calculating the net force between N particles of a k-dimensional system.  

Params (each class)
------
distance_matrix : is an instance of a DistanceMatrix.  As of now there are two implentations, either the standard distance, or a PeriodicDistance.

masses : is an array of masses for each particle

IN THE FUTURE maybe add velocity_distance_matrix or velocity_array... """
from numpy import array,vstack,sum,nan,eye
from DistanceMatrix import DistanceMatrix,PeroidicDistanceMatrix
class LeonardJonesForce(object):
  """ optional params
  ---------------
  sigma=1 : is a parameter in the Leonard Jones force

  eps=1 : is a parameter related to characteristic distance in the Leonard Jones force"""
  def __init__(self,distance_matrix,masses,sigma=1,eps=1):
    self.sigma = sigma
    self.eps = eps
    self.distance_matrix = distance_matrix
    self.masses = masses
    self.dims = len(self.masses)
  
  # this has the right signature for a Verlet integrator
  def __call__(self,x,v,t):
    return self.get_accelerations(x)

  def get_accelerations(self,x):
    distance_matrix = self.distance_matrix
    x = x.T  # Me and DistanceMatrices are fat, while containers are long
    sigma,eps,m = self.sigma,self.eps,self.masses
    dx = distance_matrix(x)		    
    r = sum(dx**2,axis=0)**.5
    r[eye(len(r),dtype='bool')] = nan	    # this is to avoid division by zero
    fmatrix = 24*eps/r*(2*(sigma/r)**12 - (sigma/r)**6)*dx
    fmatrix[:,eye(len(r),dtype='bool')] = 0 # set them back to zero
    ans = sum(fmatrix,axis=1)/m # See comment above
    return ans.T


#### Unit Testing ####
if __name__ == '__main__':
  print "--- Leonard-Jones Force ---"
  print "Periodic Boundary Conditions L = 5 in all dimensions"
  L = 5.
  distance_matrix = PeroidicDistanceMatrix(L)
  m = array([1,1,1,1])
  force = LeonardJonesForce(distance_matrix,m)
  x = array([1,2,3,4])
  y = array([3,2,1,4])
  z = array([2,3,1,4])
  p = array([x,y,z],dtype='float64')
  dpdt = eye(4)[:,:3]
  print "x = "
  print p
  print "dxdt = "
  print dpdt
  print "dvdt = "
  dvdt = force.get_accelerations(p.T)
  print "force((x,dxdt),0) = "
  print dvdt
