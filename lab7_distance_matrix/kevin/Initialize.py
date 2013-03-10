from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from numpy import linspace,sqrt,mod
from Container import Container
class ParticleInitialize(object):
  def __init__(self):
    pass
  def __call__(self,case):
    L = 10
    dims = 3
    c = Container(dims,L)
    dist = c.L[0] / 5.
    vel = dist /5.
    initialization = case
    if initialization == 'one':
	c.addParticle(0,dist,0,0,0,0,1)

    elif initialization == 'two':
	c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
	c.addParticle(dist,0.,0,-vel,0.,0.,1.)

    elif initialization == 'three':
	c.addParticle(0.,dist*sqrt(3)/2,0.,0.,-vel,0.,1.)
	c.addParticle(-dist,0.,0.,vel,0.,0.,1.0)
	c.addParticle(dist,0.,0,-vel,0.,0.,1.)

    elif initialization == 'four':
	c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
	c.addParticle(dist,0.,0,-vel,0.,0.,1.)
	c.addParticle(0.,dist,0,0.,-vel,0.,1.)
	c.addParticle(0.,-dist,0,0.,vel,0.,1.)
    elif initialization == 'six':
	c.addParticle(0.,dist,0.,0.,-vel,0.,1.)
	c.addParticle(0.,-dist,0.,0.,vel,0.,1.)
	c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
	c.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
	c.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
	c.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)
    elif initialization == 'eight':
	c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
	c.addParticle(dist,0.,0,-vel,0.,0.,1.)
	c.addParticle(0.,dist,0,0.,-vel,0.,1.)
	c.addParticle(0.,-dist,0,0.,vel,0.,1.)

	c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
	c.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
	c.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
	c.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)
    elif case == 'line':
      gamma = 1e-6
      for i in range(11):
	if i ==5:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11.,0, 1.-gamma,gamma,0,1.)
	else:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11.,0,1.,0.,0,1.) 

    elif case == 'square_lattice':
      N = 4             # Particles per row
      c.L[0] = 4.4
      c.L[1] = c.L[0]   # Extents determined by L[0] input
      d = 2.**(1/6.)    # Particle diameter
      x = linspace(d/2.,c.L[0]-d/2,N)
      y = linspace(d/2.,c.L[0]-d/2,N)
      for i in range(x.size):
	for j in range(y.size):
	  c.addParticle(x[i],y[j],0,0,0,0,1) 

    elif case == 'triangle_lattice':
      N = 8             # particles per row
      #c.L[1] = 8
      c.L[1] = sqrt(3) / 2. * c.L[0] -.2  # Set this based on L[0]
      d = 2.**(1/6.)        # diameter
      x =  linspace(-c.L[0]/2 + 3.*d/4.,c.L[0]/2. - 1.*d/4., N) # Unstaggered
      xs = linspace(-c.L[0]/2 + d/4.   ,c.L[0]/2. - 3.*d/4., N) # Staggered
      y =  linspace(-c.L[1]/2 + d/2.,c.L[1]/2  - d/2, N)
 
      for i in range(N):
        for j in range(N):
          if mod(i,2)==0:
            c.addParticle(x[j],y[i],0,0,0,0,1)
          else:
            c.addParticle(xs[j],y[i],0,0,0,0,1)
    else:
      raise ValueError("Not an option")

    distance_matrix = PeroidicDistanceMatrix(c.L)
    integrate = VerletIntegrator(.01)
    force = LeonardJonesForce(distance_matrix,c.masses)
    return c,distance_matrix,force,integrate
