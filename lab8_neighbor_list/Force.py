""" This module provides functions for calculating the net force between N particles of a k-dimensional system.  

Params (each class)
------
distance_matrix : is an instance of a DistanceMatrix.  As of now there are two implentations, either the standard distance, or a PeriodicDistance.

masses : is an array of masses for each particle

IN THE FUTURE maybe add velocity_distance_matrix or velocity_array... """
from numpy import array,vstack,sum,nan,eye,unravel_index,triu,prod,tril,size,repeat,inf,zeros,diag,ones,tile,isnan,dot
from pylab import find
#      from IPython.core.debugger import Tracer
#      debug_here = Tracer()
#      debug_here()
class SledForces(object):
  def __init__(self,dims,distance_matrix,xinit,n_sled,n_floor,load): # add sigma eps r_tol if needed
    self.leonardJones = LeonardJonesForce(dims,distance_matrix)
    self.distance_matrix = distance_matrix
    self.potential_energy = 0
    self.kinetic_energy = 0
    self.acceleration = 0
    self.pressure = 0
    self._masses = []
    self.dims = dims
    self.xinit = xinit
    self.L = array((inf,inf)) 
    self._pull_force = []
    self.n_sled = n_sled
    self.n_floor = n_floor
    self.load = load

  def append_pull(self,f):
    self._pull_force.append(f)
  @property
  def pull_force(self):
    return array(self._pull_force)
  def append_m(self,m):
    self._masses.append(m)
    self.leonardJones.append_m(m)
  @property
  def m(self):
    return array(self._masses)
  
  def __call__(self,x,v,t,calc_auxilary=True):
    load = self.load
    distance_matrix = self.distance_matrix 
    lj_force = self.leonardJones(x,v,t,calc_auxilary)
    self.potential_energy += self.leonardJones.potential_energy
    self.kinetic_energy += self.leonardJones.kinetic_energy
    dx = distance_matrix(x)
    sd_force = self.sled_springs(dx)
    pull_force = self.pull_sled(x,t,calc_auxilary)
    drag_force = self.drag_sled(x,v)
    load_force = self.load_sled(x,load)
    #print "pull force = {}".format(pull_force[-1])
    #print "time = {}".format(t)
    #return load_force + sd_force
    return (lj_force + sd_force + pull_force + drag_force + load_force)  # This is annoying, if we indeed do have bigger masses un comment .transpose()/self.m).transpose() 

  def load_sled(self,x,load):
    n_floor,n_sled = self.n_floor,self.n_sled
    force = zeros(x.shape)
    force[n_floor:,1] = -load/float(n_sled)
    return force

  def pull_sled(self,x,t,calc_auxilary=False):
    # we assume pull is on last index
    force = zeros(x.shape)
    u = x[-1,0] - self.xinit
    f = 40*(.1*t - u)
    force[-1,0] = f
    if calc_auxilary:
      self.append_pull(f)
    return force

  def drag_sled(self,x,v):
    force = zeros(x.shape)
    force[self.n_floor] = -10*v[self.n_floor]
    return force

  def sled_springs(self,dx,k=500):
    n = self.dims
    n_sled = self.n_sled
    n_floor = self.n_floor
    A = zeros((n_floor+n_sled,n_floor+n_sled,n))  # set up matrix
    B = (diag(ones(n_sled-1),k=1) + diag(ones(n_sled-2),k=2) + diag(ones(n_sled-1),k=-1) + diag(ones(n_sled-2),k=-2)) # connect the sled
    B = tile(B,(n,1,1)).transpose()  # make it n-dimensional 
    A[n_floor:,n_floor:] = B
    r = self.distance_matrix.radii(dx) 
    D = (1 - 2**(1./6.+1)/r)*dx.T*A.T*k
    D[isnan(D)] = 0
    return sum(D.T,axis=1) # calculate the force

class LeonardJonesForce(object):
  def __init__(self,dims,distance_matrix,sigma=1.,eps=1.,r_tol=1e-9):
    self.sigma = sigma
    self.eps = eps
    self.distance_matrix = distance_matrix
    self.r_tol = r_tol
    self._masses = []
    self.dims = dims
  def append_m(self,m):
    self._masses.append(m)
  def __call__(self,x,v=0,t=0,calc_auxilary=True):
    distance_matrix = self.distance_matrix
    sigma,eps,m,r_tol = self.sigma,self.eps,self._masses,self.r_tol
    dx = distance_matrix(x)		# compute distance matrix THIS IS THE MAIN TIME COMPUTATION
    r = distance_matrix.radii(dx)
    bad_diagonal = isnan(r)
    K2 = (sigma/r)**6		      # Save these constants for calculating potential energy
    K1 = K2**2 
    fmatrix = 24*eps/r*(2*K1 - K2)*dx.transpose(2,0,1)   
    K1[bad_diagonal] = 0	      # Set diagonal back to zero.  
    K2[bad_diagonal] = 0	      # Set diagonal back to zero.  
    fmatrix[:,bad_diagonal] = 0	      # Set diagonal back to zero.  
    if calc_auxilary:
      self.calculate_auxilary_forces(dx,K1,K2)
    ans = sum(fmatrix,axis=1)/m		      # Divide by masses to get acceleration
    if (abs(r) < r_tol).any():		# This is some exception handling for forces that get out of control
      self.raise_close_particle_exception(self)
    return ans.T			       

  def raise_close_particle_exception(self):
    rtol = self.rtol
    idx = find(abs(r) < r_tol)
    err_string = "Two points within 1e-8 of eachother r = \n"
    err_string += str(r) + "\n"
    err_string += "\n r{} = \n".format(array(unravel_index(idx,r.shape)))
    for x in r.ravel()[idx]:
      err_string += "{:.12f} \n".format(float(x))
    raise RuntimeError(err_string)

  def calculate_auxilary_forces(self,dx):
    self.potential_energy = sum(triu(4*eps*(K1 - K2))) # Calculate potential energy
    self.kinetic_energy = 3/2*sum(self._masses*sum(v**2,axis=1),axis=0)/2.
    N = dx.shape[0]
    d = dx.shape[2] 	     
    L = self.L[abs(self.L)>0]
    self.pressure = 0
    if N > 1:
      press_matrix = sum(fmatrix*dx.transpose(2,0,1),axis=0) # this is the "dot product"
      self.pressure = (sum(triu(press_matrix)) + N/(N-1)*self.kinetic_energy*2.)/d/prod(L)
    return self

  @property
  def L(self):
    if size(self.distance_matrix.L) == 1:
      return repeat(self.distance_matrix.L,self.dims)
    else:
      return self.distance_matrix.L
  def updateL(self,newL):
    return self.distance_matrix.updateL(newL)

#### Unit Testing ####
if __name__ == '__main__':
  from DistanceMatrix import PeroidicDistanceMatrix
  print "--- Leonard-Jones Force ---"
  print "Periodic Boundary Conditions L = 5 in all dimensions"
  L = 5.
  distance_matrix = PeroidicDistanceMatrix(L)
  m = array([1,1,1,1])
  force = LeonardJonesForce(3,distance_matrix,m)
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
  dvdt = force(p)
  print "force((x,dxdt),0) = "
  print dvdt
