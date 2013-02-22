import ODE_integrator
from numpy import array
from IPython.core.debugger import Tracer
debug_here = Tracer()
class Body(object):
  def __init__(self,initialState,mass,**kwargs):
    self.initialState = initialState
    self.mass = mass
    self.t = []
    self.pos = array(initialState[0:2])
    self.vel = array(initialState[2:])
    self.x,self.y = self.pos[0],self.pos[1]
    
  def plot(self):
    pass

class System(object):
  def __init__(self,bodies=[]):
    self.bodies = bodies

  def add_body(self,body):
    self.bodies.append(body)
    
class Force(object):
  def __init__(self):
    pass

  def __call__(self):
    pass

import csv
class InitialConditions(object):
  def __init__(self):
    pass

  def read_file(self,filename):
    with open(filename, 'r') as f:
      reader = csv.reader(f,delimiter=" ")
      numBodies = int(reader.next()[0])
      bod = []
      for row in reader:
	temp = [float(i) for i in row]  
	bod.append(Body(temp[0:-1],temp[-1]))
      return System(bodies=bod) 
  
  
