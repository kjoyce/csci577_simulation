from numpy import empty,array,zeros,vstack,hstack,inf,size,repeat,sum
from pylab import Circle,figure,gca,show
from Integrator import VerletIntegrator
#from IPython.core.debugger import Tracer
#debug_here = Tracer()

class Container(object):
  """A Container is class for computing particles simulation. An instance
    Parameters
    ----------
      dims : is the number of dimensions
      L : If the container has periodic boundary conditions, then L
      is a numpy array of lengths. """
  def __init__(self,integrator,n_sled,n_floor,dtype='float64'):
    self.dtype = dtype   # if things get too slow change this
    self._x = []  # append to this, then arrayarize.  
    self._m = []  # The @property bit makes them private
    self._v = []  # This idea came from the group Surt
    self._avg_velocities = []
    self.t = 0
    self.dims = integrator.force.dims
    self.hot_idx = None  # Hack to color the hot one
    self.integrator = integrator
    self.drag_corner_idx = n_floor
    self.pull_corner_idx = n_floor + n_sled

  @property
  def avg_sled_velocity(self):
    v = self.v[self.drag_corner_idx:]
    avg_vel = sum((v[:,0]**2 + v[:,1]**2)**.5)
    #print "avg velocity = {}".format(avg_vel)
    return avg_vel
    
  def __repr__(self):
    return_string = super(Container,self).__repr__() + "\n"
    return_string += " L: " + str(self.L) + "\n"
    for (title,e) in (("x = \n",self.x),\
		      ("v = \n",self.v),\
		      ("m = \n",self.masses)):
      return_string += title+e.__str__()+"\n"
    return return_string

  @property
  def x(self):
    #self._x = self._x % self.L # update the coordinates with each call
    return array(self._x,dtype=self.dtype) % self.L  # Note the mod L

  @property
  def v(self):
    return array(self._v,dtype=self.dtype)

  @property
  def masses(self):
    return array(self._m,dtype=self.dtype)

  @property
  def L(self):
    return self.integrator.L

  @property
  def avg_velocities(self):
    return array(self._avg_velocities)

  @property
  def pull_force(self):
    return self.integrator.force.pull_force

  def updateL(self,newL):
    return self.integrator.updateL(newL)
  def append_m(self,m):
    return self.integrator.force.append_m(m)

  def addParticle(self,*args):
    """ Adds a particle to the container and is overloaded in two
    ways.  If the particle is given as a numpy array, then there are
    two arguments, the positions and masses.  Otherwise it is a
    tuple where the last index is the mass.  Here are two examples
    that add the same particle:

    c = Container(3)
    x = array([1,2,3])
    v = zeros(3)
    m = 3.
    c.addParticle(x,v,m)

    c.addParticle(1,2,3,0,0,0,3.)"""
    if len(args) == 3: 
      self._x.append(args[0])
      self._v.append(args[1])
      self.append_m(args[2])
    elif len(args)%2 and len(args)>1: # This must be odd
      dim = len(args)/2
      self._x.append(args[0:dim])
      self._v.append(args[dim:-1])
      self.append_m(args[-1])
    else:
      raise ValueError('must pass odd numer of args greater than 2')

  def integrate(self):
    self.t += self.integrator.dt
    (dx,dv) = self.integrator(self.x,self.v,self.t)
    move_from = self.drag_corner_idx
    self._x[move_from:] = (self.x + dx)[move_from:]
    self._v[move_from:] = (self.v + dv)[move_from:]
    self._avg_velocities.append(self.avg_sled_velocity)

  def etargetni(self):
    (dx,dv) = self.integrator.backward(self.x,self.v,self.t)
    self._x = self.x + dx
    self._v = self.v + dv
    
if __name__ == '__main__':
  L = array([2,3,4])
  c = Container(3,L)
  print "L = {}".format(L)
  x = array([1,3,4])
  v = zeros(3)
  m = 3.
  print "adding x = {}...".format(x)
  c.addParticle(x,v,m)
  print c
  args = range(7)
  print "adding with tuple {}".format(args)
  c.addParticle(*args)
  print c
  
