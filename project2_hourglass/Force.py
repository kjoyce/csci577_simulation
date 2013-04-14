from numpy import array,vstack,sum,nan,eye,unravel_index,triu,prod,tril,size,repeat,inf,zeros,diag,ones,tile,isnan,dot,pi,sin,cos
from pylab import find
from scipy.sparse import csc_matrix
#      from IPython.core.debugger import Tracer
#      debug_here = Tracer()
#      debug_here()
class HourGlassForce(object):
  def __init__(self,distance_matrix,angle=pi/3,holewidth=2):
    self.distance_matrix = distance_matrix
    self.dims = 2
  
  def __call__(self,x,v,t,sparse_signature=None):
    if sparse_signature is None:
      dx,dy   = self.dense_dmatrix(x)
      dvx,dvy = self.dense_dmatrix(v)
      dr = (dx**2 + dy**2)**.5
      sparse_structure = (dr < 2.**(1./6.))*(dr > 0)
      sparse_signature = csc_matrix(sparse_structure,dtype='float64')
      self.sparse_signature = sparse_signature
      sparcify = lambda dx: csc_matrix((dx[sparse_signature.nonzero()],sparse_signature.indices,sparse_signature.indptr))
      (dx,dy,dr,dvx,dvy) =  map(sparcify, (dx,dy,dr,dvx,dvy))
    else:
      update = lambda x: self.sparse_dmatrix(x,sparse_signature)
      (dx,dy,dvx,dvy) = map(update, (x.T[0],x.T[1],v.T[0],v.T[1]))
      dr = dx.multiply(dx) + dy.multiply(dy)
      dr.data = dr.data**.5

    lj_force = self.lennard_jones_force(dx,dy,dr)
    ds_force = self.dissipation_force(dx,dy,dr,dvx,dvy)
    gravity = -2.*ones(lj_force.shape)
    gravity[0] = 0
    total_force = (lj_force + ds_force + gravity)
    print "d_force: {}".format(ds_force)
    return total_force.T

  @classmethod
  def dense_dmatrix(cls,p):
    dp = tile(p,(len(p),1,1))
    return (dp.transpose((1,0,2)) - dp).T

  @classmethod
  def sparse_dmatrix(cls,x,sparse_signature):
    sp_dx = cls.sparse_tile(x,sparse_signature)
    return (sp_dx.transpose() - sp_dx)

  @classmethod
  def sparse_tile(cls,x_new,sp_dx):
    new_sp_dx = sp_dx.astype('float64')    # make a copy, just for illustrative purposes
    for i,j in zip(sp_dx.indptr[:-1],sp_dx.indptr[1:]):
      new_sp_dx.data[i:j] = x_new[sp_dx.indices[i:j]]
    return new_sp_dx

  def lennard_jones_force(self,dx,dy,dr):
    fxmatrix = dr.copy()  # we assume that dx,dy, and r have same sparcity
    fymatrix = dr.copy()
    K1 = dr.data**(-6.)
    K2 = K1**2
    fxmatrix.data = 24./dr.data*(2*K1 - K2)
    fxmatrix = fxmatrix.multiply(dx)
    fymatrix.data = 24./dr.data*(2*K1 - K2)
    fymatrix = fymatrix.multiply(dy)
    return array(vstack([fxmatrix.sum(0),fymatrix.sum(0)]))

  def dissipation_force(self,dx,dy,dr,dvx,dvy,gamma=10):
    fxmatrix = dr.copy()  
    fxmatrix.data = fxmatrix.data**(-2.)
    fymatrix = dr.copy()
    fymatrix.data = fymatrix.data**(-2.)
    dot_prod = dvx.multiply(dx) + dvy.multiply(dy)
    fxmatrix = fxmatrix.multiply(dot_prod).multiply(dx)
    fymatrix = fymatrix.multiply(dot_prod).multiply(dy)
    return -gamma*array(vstack([fxmatrix.sum(0),fymatrix.sum(0)]))

  def hourglass_force(self,x,v,gamma=20):
    pass

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
  from DistanceMatrix import DistanceMatrix
  from pylab import *
  distance_matrix = DistanceMatrix()
  hg_force = HourGlassForce(2,distance_matrix)
  sig = .5*2.**(1./6.)
  x = sig*array([[0,0],[1,0],[5,0],[6,0],[7,0],[1.5,(sqrt(3)+1)/2]])
  v = zeros(x.shape)
  f = hg_force(x,v,0)
  print "force: {} \n sp_sig: {}".format(f,hg_force.sparse_signature)
  f = hg_force(x,v,1,hg_force.sparse_signature)
  print "force: {} \n sp_sig: {}".format(f,hg_force.sparse_signature)

