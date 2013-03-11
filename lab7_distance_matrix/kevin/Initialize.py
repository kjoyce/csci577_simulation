from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from numpy import linspace,sqrt,mod
from pylab import rand,randint
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
    xlim = (-5.,15.)
    ylim = (-5.,15.)
    pot_energy_lim = (-10,30)
    kin_energy_lim = (40,70)
    tot_energy_lim = (40,70)
    if initialization == 'one':
      c.addParticle(0,dist,0,0,0,0,1)

    elif initialization == 'two':
      c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
      c.addParticle(dist,0.,0,-vel,0.,0.,1.)
      pot_energy_lim = (-1,1)
      kin_energy_lim = (-1,2)
      tot_energy_lim = (0,2)
      pressure_lim   = (-1,1)

    elif initialization == 'three':
      c.addParticle(0.,dist*sqrt(3)/2,0.,0.,-vel,0.,1.)
      c.addParticle(-dist,0.,0.,vel,0.,0.,1.0)
      c.addParticle(dist,0.,0,-vel,0.,0.,1.)
      pot_energy_lim = (-2,1)
      kin_energy_lim = (-1,3)
      tot_energy_lim = (0,2)
      pressure_lim   = (-1,1)

    elif initialization == 'four':
      c.addParticle(-dist,0.,0.,vel,0.,0.,1.)
      c.addParticle(dist,0.,0,-vel,0.,0.,1.)
      c.addParticle(0.,dist,0,0.,-vel,0.,1.)
      c.addParticle(0.,-dist,0,0.,vel,0.,1.)
      pot_energy_lim = (-5,5)
      kin_energy_lim = (-.5,10)
      tot_energy_lim = (0,4)
      pressure_lim   = (-4,4)

    elif initialization == 'six':
      energy_lim = 80
      c.addParticle(0.,dist,0.,0.,-vel,0.,1.)
      c.addParticle(0.,-dist,0.,0.,vel,0.,1.)
      c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
      c.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
      c.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
      c.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)
	
    elif initialization == 'eight':
      energy_lim = 80
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
      pot_energy_lim = (-10,30)
      kin_energy_lim = (40,70)
      xlim = (-5,15)
      ylim = xlim
      for i in range(11):
	if i ==5:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11.,0, 1.-gamma,gamma,0,1.)
	else:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11.,0,1.,0.,0,1.) 

    elif case == 'square_lattice':
      N = 4             # Particles per row
      energy_lim = 30
      xlim = (-2,6)
      ylim = xlim
      c.L[0] = 4.4
      c.L[1] = c.L[0]   # Extents determined by L[0] input
      d = 2.**(1/6.)    # Particle diameter
      x = linspace(d/2.,c.L[0]-d/2,N)
      y = linspace(d/2.,c.L[0]-d/2,N)
      for i in range(x.size):
	for j in range(y.size):
	  c.addParticle(x[i],y[j],0,0,0,0,1) 
    
    elif case == 'hot_square_lattice':
      N = 4             # Particles per row
      energy_lim = 400
      xlim = (-2,6)
      ylim = xlim
      c.L[0] = 4.4
      c.L[1] = c.L[0]   # Extents determined by L[0] input
      d = 2.**(1/6.)    # Particle diameter
      x = linspace(d/2.,c.L[0]-d/2,N)
      y = linspace(d/2.,c.L[0]-d/2,N)
      k = 0
      hot_idx = randint(N**2)
      for i in range(x.size):
	for j in range(y.size):
	  vx = (hot_idx == i*j)*rand()
	  vy = (hot_idx == i*j)*rand()
	  c.addParticle(x[i],y[j],0,vx,vy,0,1) 
      
    elif case == 'triangle_lattice':
      ylim = (-1,9)
      xlim = (-1,10)
      N = 8             # particles per row
      energy_lim = 400
      c.L[1] = 8
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
    return c,distance_matrix,force,integrate,xlim,ylim,pot_energy_lim,kin_energy_lim,tot_energy_lim,pressure_lim
