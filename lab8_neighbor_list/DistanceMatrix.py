from numpy import array,tile,size,inf
#from IPython.core.debugger import Tracer
#debug_here = Tracer()
class DistanceMatrix(object):
  def __init__(self):
    self._L = inf
  def __call__(self,p):
    dp = self.all_pair_matrix(p)
    self.dp = self.signed_distance(dp,dp.transpose((1,0,2))) # tranposes dimensions according to 0->1 1->0 2->2
    return self.dp
  def all_pair_matrix(self,p):
    return tile(p,(len(p),1,1))  
  def signed_distance(self,x,y):
    return x - y
  @property
  def L(self):
    return self._L
  def updateL(self,newL):
    self._L = newL

class PeroidicDistanceMatrix(DistanceMatrix):
  def __init__(self,L):
    self.updateL(L)
  # Think hard about why this works.  This is the map (x,y) -> d
  # where d is the one dimensional toroidal distance. It is
  # essentialy one phase (y dimension) of a ''flat'' sine wave (x
  # dimension).  Remember, a%b = a mod b 
  def signed_distance(self,x,y):
    L = self.L
    flat_wave = lambda x,L: (x-.5*L)%L-.5*L
    return flat_wave(x%L-y%L,L)

#### Unit Testing ####
if __name__ == '__main__':
  from numpy import size,vstack,mgrid,zeros,unravel_index
  from matplotlib.pyplot import *
  print '--Testing--'
  x = array([1,2,3,4])
  y = array([3,2,1,4])
  z = array([2,3,1,4])
  p = vstack([x,y,z])
  distance_matrix = DistanceMatrix()
  dp = distance_matrix(p.T)
  (dx,dy,dz) = dp.T
  print "x = {}".format(x)
  print "y = {}".format(y)
  print "z = {}".format(z)
  print "dx = \n{}".format(dx)
  print "dy = \n{}".format(dy)
  print "dz = \n{}".format(dz)

  print '--Now for something crazy--'
  from itertools import permutations
  x = array([x for x in permutations(range(3))])
  print "6 dimesional array \nx = \n{}".format(x)
  dx = distance_matrix(x.T)
  print "dx =\n {}".format(dx)

  print '--Period Boundary Conditions--'
  L = 1.
  num_pts = 20.
  pt_idx = (int(.15*num_pts),int(.30*num_pts))

  delta = num_pts**-1
  idx = pt_idx[0]*num_pts + pt_idx[1]  # C style raveling
  pdistance_matrix = PeroidicDistanceMatrix(L)
  x,y = mgrid[0:L:delta*L,0:L:delta*L]
  p = vstack((x.ravel(),y.ravel()))  
  dp = pdistance_matrix(p.T)
  picture = zeros((num_pts,num_pts))
  picture[unravel_index(range(int(num_pts**2)),(num_pts,num_pts))] = dp[:,idx,0]**2 + dp[:,idx,1]**2
  #plot( pt_idx[1], pt_idx[0] , 'ob', markersize=10)  # this is if you want fancy interpolation
  #imshow(picture,origin='lower',interpolation="None")
  plot( x[pt_idx], y[pt_idx] , 'ob', markersize=10)
  contourf(x,y,picture)
  colorbar()
  title("Domain: $[0.0,{}]$, Distance from $({},{})$".format(L,x[pt_idx],y[pt_idx]))

  show()

  L = 5.
  print "L = {}".format(L)
  x = array([1,2,3,4])
  y = array([3,2,1,4])
  z = array([2,3,1,4])
  p = vstack([x,y,z])
  pdistance_matrix = PeroidicDistanceMatrix(L)
  dp = pdistance_matrix(p.T)
  (dx,dy,dz) = dp.T
  print "x = {}".format(x)
  print "y = {}".format(y)
  print "z = {}".format(z)
  print "dx = \n{}".format(dx)
  print "dy = \n{}".format(dy)
  print "dz = \n{}".format(dz)
