from numpy import empty,array,zeros,vstack,hstack,inf,size,repeat
from pylab import Circle,figure,gca,show
#from IPython.core.debugger import Tracer
#debug_here = Tracer()

class Container(object):
  """A Container is class for computing particles simulation. An instance
    Parameters
    ----------
      dims : is the number of dimensions
      L : If the container has periodic boundary conditions, then L
      is a numpy array of lengths. """
  def __init__(self,dims,L=inf,dtype='float64'):
    self.dtype = dtype   # if things get too slow change this
    if size(L) == 1:
      self.L = repeat(L,dims)
    else:
      self.L = L
    self._x = []  # append to this, then arrayarize.  
    self._m = []  # The @property bit makes them private
    self._v = []  # This idea came from the group Surt
    self.dims = dims
    self.circles = []

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
    return array(self._x,dtype=self.dtype) % self.L  # Note the mod L

  @property
  def v(self):
    return array(self._v,dtype=self.dtype)

  @property
  def masses(self):
    return array(self._m,dtype=self.dtype)

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
      self._m.append(args[2])
    elif len(args)%2 and len(args)>1: # This must be odd
      dim = len(args)/2
      self._x.append(args[0:dim])
      self._v.append(args[dim:-1])
      self._m.append(args[-1])
    else:
      raise ValueError('must pass odd numer of args greater than 2')

  def __make_circle(self, xy, radius, ax, color="lightsteelblue", facecolor="green", alpha=.6):
    """ add a circle to ax= or current axes
    """
    e = Circle( xy=xy[::-1], radius=radius ) # notice the stupid CS reversed coordinates
    ax.add_artist(e)
    e.set_clip_box(ax.bbox)
    e.set_edgecolor( color )
    e.set_linewidth(3)
    e.set_facecolor( facecolor )  # "none" not None
    e.set_alpha( alpha )
    e.set_animated( True )
    print "HERE: "
    print e
    return e

  def __set_all_circles(self,ax):
    for x in self.x:
      self.circles.append(self.__make_circle(x,.1,ax))

  def __update_all_circles(self,ax):
    for i in range(len(self.circles)):
      self.circles[i].update_from(self.__make_circle(self.x[i],.1,ax))

  def draw(self,ax):
    ax.set_xlim((0,self.L[0]))
    ax.set_ylim((0,self.L[1]))
    if not(self.circles):
      self.__set_all_circles(ax)
    else:
      self.__update_all_circles(ax)
      
  def integrate(self,dx,dv):
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
  c = Container(2,5)
  c.addParticle(2,2,0,0,1.)
  c.addParticle(3,6,0,0,1.)
  print c
  figure(figsize=(c.L[0]*1.3,c.L[1]*1.3))
  ax = gca()
  c.draw(ax)
  ax.set_aspect('equal')
  show()
  
