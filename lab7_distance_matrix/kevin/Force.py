""" This module provides functions for calculating the net force between N particles of a k-dimensional system.  

Params (each class)
------
distance_matrix : is an instance of a DistanceMatrix.  As of now there are two implentations, either the standard distance, or a PeriodicDistance.

masses : is an array of masses for each particle

IN THE FUTURE maybe add velocity_distance_matrix or velocity_array... """
from numpy import array,vstack,sum,nan,eye,unravel_index,triu,prod,tril
from pylab import find
from DistanceMatrix import DistanceMatrix,PeroidicDistanceMatrix
class LeonardJonesForce(object):
  """ optional params
  ---------------
  sigma=1 : is a parameter in the Leonard Jones force

  eps=1 : is a parameter related to characteristic distance in the Leonard Jones force"""
  def __init__(self,distance_matrix,masses,sigma=1.,eps=1.,r_tol=1e-9):
    self.sigma = sigma
    self.eps = eps
    self.distance_matrix = distance_matrix
    self.masses = masses
    self.dims = len(self.masses)
    self.r_tol = r_tol
  
  # this has the right signature for a Verlet integrator
  def __call__(self,x,v,t,calc_auxilary=True):
    distance_matrix = self.distance_matrix
    sigma,eps,m,r_tol = self.sigma,self.eps,self.masses,self.r_tol
    dx = distance_matrix(x)		# compute distance matrix THIS IS THE MAIN TIME COMPUTATION
    r = sum(dx**2,axis=2)**.5		# precalculate radii
    r[eye(len(r),dtype='bool')] = nan	# this is to avoid division by zero
#    r[ r==0 ] = 1			# This is some exception handling for forces that get out of control
    if (abs(r) < r_tol).any():
      idx = find(abs(r) < r_tol)
      err_string = "Two points within 1e-8 of eachother r = \n"
      err_string += str(r) + "\n"
      err_string += "\n r{} = \n".format(array(unravel_index(idx,r.shape)))
      for x in r.ravel()[idx]:
	err_string += "{:.12f} \n".format(float(x))
      raise ValueError(err_string)
    bad_diagonal = eye(len(r),dtype='bool')
    K1 = (sigma/r)**12 
    K2 = (sigma/r)**6			      # Save these constant for calculating potential energy
    fmatrix = 24*eps/r*(2*K1 - K2)*dx.transpose(2,0,1)  # 
    fmatrix[:,bad_diagonal] = 0		      # Set diagonal back to zero
    K1[bad_diagonal] = 0		      # Set diagonal back to zero
    K2[bad_diagonal] = 0		      # Set diagonal back to zero
    if calc_auxilary:
#      from IPython.core.debugger import Tracer
#      debug_here = Tracer()
#      debug_here()
      self.potential_energy = sum(triu(4*eps*(K1 - K2))) # Calculate potential energy
      self.kinetic_energy = 3/2*sum(self.masses*sum(v**2,axis=1),axis=0)/2.
      N = dx.shape[0]
      d = dx.shape[2] -1	      ######## HACKED
      L = self.distance_matrix.L[0:1] ######## FOR FORCED 2D
      self.pressure = 0
      if N > 1:
	press_matrix = sum(fmatrix*dx.transpose(2,0,1),axis=0) # this is the "dot product"
	self.pressure = (sum(triu(press_matrix)) + N/(N-1)*self.kinetic_energy*2.)/d/prod(L)
    ans = sum(fmatrix,axis=1)/m		      # Divide by masses to get acceleration
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
  p = array([x,y,z],dtype='float64').T
  dpdt = eye(4)[:,:3]
  print "x = "
  print p
  print "dxdt = "
  print dpdt
  print "dvdt = "
  dvdt = force.get_accelerations(p)
  print "force((x,dxdt),0) = "
  print dvdt
