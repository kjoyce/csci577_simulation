from numpy import empty,array,zeros,hstack
class Container(object):
  """A Container is class for computing particles simulation. An instance
  requires the shape of the vectors to process. For a periodic domain include
  the keyword argument L.  

    Parameters
    ----------
      shape : a tuple of integers that describes the dimension of
      the spacial variable.  

      L : If the container has periodic boundary conditions, then L
      is a numpy array of lengths.  It must be true that L.shape =
      shape, otherwise throw some an exception. """
  def __init__(self,shape,**kwargs):
    dattype = 'float64'   # if things get too slow change this
    self.x,self.v,self.masses = [empty(0,dtype=dattype) for i in range(3)]
    if 'L' in kwargs:
      self.L = L
      self.x = empty(L.shape)
      if kwargs['L'].shape != shape:
	raise ValueError('arrays must have the same number of dimensions')


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

    c.addParticle(1,2,3,0,0,0,1,3.)"""
    if len(args) == 3: 
      x,v,m = self.x,self.v,self.masses
      self.x,self.v,self.masses = [hstack([e,args[i]]) for (e,i) in zip((x,v,m),range(3))] 
    elif len(args) % 2: # This must be odd
      dim = len(args)/2
    else:
      raise ValueError('must pass odd numer of args')
      
  def __repr__(self):
    return_string = ""
    for (title,e) in (("x = ",self.x),\
		      ("v = ",self.v),\
		      ("m = ",self.masses)):
      return_string += title+e.__str__()+"\n"
    return return_string

if __name__ == '__main__':
  c = Container(3)
  x = array([1,2,3])
  v = zeros(3)
  m = 3.
  c.addParticle(x,v,m)
  print c
  c.addParticle(1,2,3,0,0,0,3.)
  print c
