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
  
  def __call__(self,x,v,t,calc_auxilary=True,sparse_signature=None):
    if sparse_signature is None:
      dx,dy   = self.dense_dmatrix(x)
      dvx,dvy = self.dense_dmatrix(v)
      dr = (dx**2 + dy**2)**.5
      sparse_structure = (dr < 2.**(1./6.))*(dr > 0)
      sparse_signature = csc_matrix(sparse_structure)
      sparcify = lambda dx: csc_matrix((dx[sparse_signature.nonzero()],sparse_signature.indices,sparse_signature.indptr))
      (dx,dy,dr,dvx,dvy) =  map(sparcify, (dx,dy,dr,dvx,dvy))
    else:
      update = lambda x: self.sparse_tile(x,sparse_signature)
      (dx,dy,dvx,dvy) = map(update, (x.T[0],x.T[1],v.T[0],v.T[1]))
      dr = (dx**2 + dy**2)**.5

    lj_force = self.lennard_jones_force(dx,dy,dr)
    ds_force = self.dissipation_force(dx,dy,dr,dvx,dvy)
    return (lj_force).T,sparse_signature

  @classmethod
  def dense_dmatrix(cls,p):
    dp = tile(x,(len(x),1,1))
    return (dp.transpose((1,0,2)) - dp).T

  @classmethod
  def sparse_tile(cls,x_new,sp_dx):
    new_sp_dx = sp_dx.copy()    # make a copy, just for illustrative purposes
    for i,j in zip(sp_dx.indptr[:-1],sp_dx.indptr[1:]):
      new_sp_dx.data[i:j] = x_new[sp_dx.indices[i:j]]
    return new_sp_dx

  def lennard_jones_force(self,dx,dy,dr):
    fxmatrix = dr.copy()  # we assume that dx,dy, and r have same sparcity
    fymatrix = dr.copy()
    K1 = dr.data**(-6.)
    K2 = K1**2
    fxmatrix.data = 24./dr.data*(2*K1 - K2)*dx.data   
    fymatrix.data = 24./dr.data*(2*K1 - K2)*dy.data   
    return array(vstack([fxmatrix.sum(0),fymatrix.sum(0)]))

  def dissipation_force(self,dx,dy,dr,dvx,dvy,gamma=20):
    fxmatrix = dr.copy()  # we assume that dx,dy, and r have same sparcity
    fymatrix = dr.copy()
    dot_prod = (dvx.data*dx.data + dvy.data*dx.data)
    fxmatrix.data = -gamma*dot_prod*dx.data/dr.data**2
    fymatrix.data = -gamma*dot_prod*dy.data/dr.data**2
    return array(vstack([fxmatrix.sum(0),fymatrix.sum(0)]))

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
  sig = 2.**(1./6.)
  x = array([[0,0],[1,0],[5,0],[6,0],[7,0],[1.5,(sqrt(3)+1)/2]])
  v = zeros(x.shape)
  f,sp_sig = hg_force(x,v,0)
  print "force: {} \n sp_sig: {}".format(f,sp_sig)
  f,sp_sig = hg_force(x,v,1,sp_sig)
  print "force: {} \n sp_sig: {}".format(f,sp_sig)

