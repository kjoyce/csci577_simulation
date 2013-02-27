import csv
from scipy.integrate import odeint  # for integrate.odeint
from numpy import *
from pylab import *
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
  def __init__(self,bodies):
    self.bodies = bodies
    n = size(bodies)
    pos = []
    vel = []
    for i in range(n):
      pos.append(bodies[i].pos)
      vel.append(bodies[i].vel)
    self.init = array([pos,vel]).T.flatten(1)

class Force(object):
  def __init__(self,system,constant):
    self.constant = constant
    self.system = system
    self.nbodies = size(system.bodies)
  
  def __call__(self,x,t):
    pairs = zeros((2,self.nbodies,self.nbodies))
    bod = self.system.bodies
    G = self.constant
    for i in range(size(bod)):
      x1 = x[2*i:2*i+2]
      m1 = bod[i].mass
      for j in range(i+1,size(bod)):
	x2 = x[2*j:2*j+2]
	m2 = bod[j].mass
	pairs[:,i,j] = self.interaction(x1,m1,x2,m2)
    accelerationsx = sum(pairs[0] - pairs[0].T,axis=1)
    accelerationsy = sum(pairs[1] - pairs[1].T,axis=1)
    accelerations = array([accelerationsx,accelerationsy]).flatten(1)
    return array([x[2*size(bod):],accelerations]).T.flatten(1)
    
  def interaction(self,x1,m1,x2,m2):
    G = self.constant
    r=norm(x2-x1)
    xdist=(x1[0]-x2[0])
    ydist=(x1[1]-x2[1])
    ax=-G*m2*xdist/r**3
    ay=-G*m2*ydist/r**3
    return array([ax,ay])



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
  
G =1
files = ('euler_init.txt','montgomery_init.txt','lagrange_init.txt')
for fname in files:
  reader = InitialConditions()
  sys_1 = reader.read_file(fname)
  force = Force(sys_1,G)
  t = linspace(0.,100.,1000)
  x = odeint(force,sys_1.init,t)
  figure()
  plot(x[:,0],x[:,1])
  plot(x[:,2],x[:,3])
  plot(x[:,4],x[:,5])
  show()
