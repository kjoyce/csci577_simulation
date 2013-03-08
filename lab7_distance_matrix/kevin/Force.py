""" This module provides functions for calculating the net force between N particles of a k-dimensional system.  

Params (each class)
------
distance_matrix : is an instance of a DistanceMatrix.  As of now there are two implentations, either the standard distance, or a PeriodicDistance.

masses : is an array of masses for each particle

IN THE FUTURE maybe add velocity_distance_matrix or velocity_array... """
from numpy import array,vstack,sum,nan,eye,unravel_index
from pylab import find
from DistanceMatrix import DistanceMatrix,PeroidicDistanceMatrix
class LeonardJonesForce(object):
  """ optional params
  ---------------
  sigma=1 : is a parameter in the Leonard Jones force

  eps=1 : is a parameter related to characteristic distance in the Leonard Jones force"""
  def __init__(self,distance_matrix,masses,sigma=1.,eps=1.):
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
    x = x.T				    # DistanceMatrices are fat, while (which this returns to) containers are long and slim
    sigma,eps,m = self.sigma,self.eps,self.masses
    dx = distance_matrix(x)		    # compute distance matrix THIS IS THE MAIN TIME COMPUTATION
    r = sum(dx**2,axis=0)**.5		    # precalculate radii
    r[eye(len(r),dtype='bool')] = nan	    # this is to avoid division by zero
#    r[ r==0 ] = 1							    This is some exception handling for forces that get out of control
    if (abs(r) < 1e-8).any():
      idx = find(abs(r) < 1e-8)
      err_string = "Two points within 1e-8 of eachother r = \n"
      err_string += str(r) + "\n"
      err_string += "\n r{} = \n".format(array(unravel_index(idx,r.shape)))
      for x in r.ravel()[idx]:
	err_string += "{:.12f} \n".format(float(x))
      raise ValueError(err_string)
    K = (2*(sigma/r)**12 - (sigma/r)**6)      # save this constant for calculation potential energy
    fmatrix = 24*eps/r*K*dx		      # calculate force
    fmatrix[:,eye(len(r),dtype='bool')] = 0   # set them back to zero
    self.potential_energy = sum(4*eps*K[K>0]) # calculate potential energy
    ans = sum(fmatrix,axis=1)/m		      # Divide by masses to get acceleration
    return ans.T			      # See fat/slim comment

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
