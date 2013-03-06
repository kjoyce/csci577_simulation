from numpy import linspace
from Container import Container
class ParticleInitialize(object):
  def __init__(self):
    pass
  def __call__(self,case,c):
    if case == 'VerletTest':
      pass

    elif case == 'line':
      gamma = 1e-6
      for i in range(11):
	if i ==5:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11., 1.-gamma,gamma,1.)
	else:
	  c.addParticle(c.L[0] / 2., (i-.5) * c.L[1] / 11.,1.,0.,1.) 

    elif case == 'square_lattice':
      N = 8             # Particles per row
      c.L[1] = c.L[0]       # Extents determined by L[0] input
      d = 2.**(1/6.)    # Particle diameter
      x = linspace(-c.L[0]/2+d/2.,c.L[0]/2-d/2,N)
      y = linspace(-c.L[0]/2+d/2.,c.L[0]/2-d/2,N)
      for i in range(x.size):
	for j in range(y.size):
	  c.addParticle(x[i],y[j],0,0,0,0,1) 

    elif case == 'triangle_lattice':
      N = 8             # particles per row
      c.L[1] = sqrt(3) / 2. * c.L[0]  # Set this based on L[0]
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
    return c
