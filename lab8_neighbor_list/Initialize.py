from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from numpy import linspace,sqrt,mod
from pylab import rand,randint
from Container import Container
class ParticleInitialize(object):
  def __init__(self,case):
    L = 10
    dims = 2
    dt = .01
    self.distance_matrix = PeroidicDistanceMatrix(L)
    self.force		 = LeonardJonesForce(dims,self.distance_matrix)
    self.integrate	 = VerletIntegrator(dt,self.force)
    self.c		 = Container(self.integrate)
    c = self.c
    dist = c.L[0] / 5.
    vel = dist /5.
    initialization = case
    xlim = (-5.,15.)
    ylim = (-5.,15.)
    pot_energy_lim = (-10,30)
    kin_energy_lim = (40,70)
    tot_energy_lim = (40,70)

    if initialization == 'one':
      c.addParticle(0,dist,0,0,1)

    elif initialization == 'two':
      c.addParticle(-dist,0.,vel,0.,1.)
      c.addParticle(dist,0.,-vel,0.,1.)
      pot_energy_lim = (-1,1)
      kin_energy_lim = (-1,2)
      tot_energy_lim = (0,2)
      pressure_lim   = (-1,1)

    elif initialization == 'three':
      c.addParticle(0.,dist*sqrt(3)/2.,0.,-vel,1.)
      c.addParticle(-dist,0.,vel,0.,1.0)
      c.addParticle(dist,0.,-vel,0.,1.)
      pot_energy_lim = (-2,1)
      kin_energy_lim = (-1,3)
      tot_energy_lim = (0,2)
      pressure_lim   = (-1,1)

    elif initialization == 'four':
      c.addParticle(-dist,0.,vel,0.,1.)
      c.addParticle(dist,0.-vel,0.,1.)
      c.addParticle(0.,dist,0.,-vel,1.)
      c.addParticle(0.,-dist,0.,vel,1.)
      pot_energy_lim = (-5,5)
      kin_energy_lim = (-.5,10)
      tot_energy_lim = (0,4)
      pressure_lim   = pot_energy_lim

    elif initialization == 'six':
      energy_lim = 80
      c.addParticle(0.,dist,0.,-vel,1.)
      c.addParticle(0.,-dist,0.,vel,1.)
      c.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
      c.addParticle(-dist/sqrt(2),dist/sqrt(2),vel/sqrt(2),-vel/sqrt(2),1.)
      c.addParticle(-dist/sqrt(2),-dist/sqrt(2),vel/sqrt(2),vel/sqrt(2),1.)
      c.addParticle(dist/sqrt(2),-dist/sqrt(2),-vel/sqrt(2),vel/sqrt(2),0.,1.)
      pot_energy_lim = (-10,10)
      kin_energy_lim = (-.5,10)
      tot_energy_lim = (-10,10)
      pressure_lim   = pot_energy_lim
	
    elif initialization == 'eight':
      energy_lim = 80
      c.addParticle(-dist,0.,vel,0.,1.)
      c.addParticle(dist,0,-vel,0.,1.)
      c.addParticle(0.,dist,0.,-vel,1.)
      c.addParticle(0.,-dist,0.,vel,1.)

      c.addParticle(dist/sqrt(2),dist/sqrt(2),-vel/sqrt(2),-vel/sqrt(2),1.)
      c.addParticle(-dist/sqrt(2),dist/sqrt(2),vel/sqrt(2),-vel/sqrt(2),1.)
      c.addParticle(-dist/sqrt(2),-dist/sqrt(2),vel/sqrt(2),vel/sqrt(2),1.)
      c.addParticle(dist/sqrt(2),-dist/sqrt(2),-vel/sqrt(2),vel/sqrt(2),1.)
      pot_energy_lim = (-10,10)
      kin_energy_lim = (-.5,10)
      tot_energy_lim = (-10,10)
      pressure_lim   = (-15,15)

    elif case == 'line':
      gamma = 1e-6
      pot_energy_lim = (-10,100)
      kin_energy_lim = (-.5,100)
      xlim = (-1,11)
      ylim = xlim
      tot_energy_lim = (0,100)
      pressure_lim   = (0,30)
      for i in range(11):
	if i ==5:
	  c.addParticle(c.L[0]/2.,(i-.5)*c.L[1]/11.,1.-gamma,gamma,1.)
	else:
	  c.addParticle(c.L[0]/2.,(i-.5)*c.L[1]/11.,1.,0.,1.) 

      
    elif case == 'square_lattice' or case == 'crunch_square_lattice' or case == 'hot_square_lattice' or case =='square_lattice8x8':
      if case == 'hot_square_lattice':
	c.hot_idx = randint(64)
      N = 8             # Particles per row
      xlim = (-1,10)
      ylim = xlim
      if case == 'crunch_square_lattice':
	pot_energy_lim = (-300,1000)
	kin_energy_lim = (-.5,1000)
	tot_energy_lim = (-300,1000)
	pressure_lim   = (-200,1000)
      else:
	pot_energy_lim = (-300,100)
	kin_energy_lim = (-.5,40)
	tot_energy_lim = (-300,100)
	pressure_lim   = (-200,100)
      c.updateL(8)
      d = 2.**(1/6.)    # Particle diameter
      x = linspace(d/2.,c.L[0]-d/2,N)
      y = linspace(d/2.,c.L[0]-d/2,N)
      for i in range(x.size):
	for j in range(y.size):
	  if i*8 + j == c.hot_idx:
	    (vx,vy) = (rand(),rand())
	    print "(vx,vy) = ({},{}) at ({},{})".format(vx,vy,i,j)
	    c.addParticle(x[i],y[j],vx,vy,1) 
	  else:
	    c.addParticle(x[i],y[j],0,0,1) 
    
    elif case == 'triangle_lattice' or case == 'crunch_triangle_lattice':
      ylim = (-1,9)
      xlim = (-1,11)
      N = 8             # particles per row
      if case == 'crunch_triangle_lattice':
	pot_energy_lim = (-300,1000)
	kin_energy_lim = (-.5,1000)
	tot_energy_lim = (-300,1000)
	pressure_lim   = (-200,1000)
      else:
	pot_energy_lim = (-300,100)
	kin_energy_lim = (-.5,15)
	tot_energy_lim = (-300,100)
	pressure_lim   = (-200,100)
      c.L[0] = 8.8
      c.L[1] = sqrt(3) / 2. * c.L[0] -.2  # Set this based on L[0]
      d = 2.**(1/6.)        # diameter
      x =  linspace(-c.L[0]/2 + 3.*d/4.,c.L[0]/2. - 1.*d/4., N) # Unstaggered
      xs = linspace(-c.L[0]/2 + d/4.   ,c.L[0]/2. - 3.*d/4., N) # Staggered
      y =  linspace(-c.L[1]/2 + d/2.,c.L[1]/2  - d/2, N)
 
      for i in range(N):
        for j in range(N):
          if mod(i,2)==0:
            c.addParticle(x[j],y[i],0,0,1)
          else:
            c.addParticle(xs[j],y[i],0,0,0,0,1)
    else:
      raise ValueError("Not an option")

    self.xlim            = xlim
    self.ylim            = ylim
    self.pot_energy_lim  = pot_energy_lim
    self.kin_energy_lim  = kin_energy_lim
    self.tot_energy_lim  = tot_energy_lim
    self.pressure_lim    = pressure_lim
  def __call__(self):
    return self.c,self.distance_matrix,self.force,self.integrate,self.xlim,self.ylim,self.pot_energy_lim,self.kin_energy_lim,self.tot_energy_lim,self.pressure_lim


