from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from Initialize import ParticleInitialize
from pylab import figure,clf,plot,title,show,gca,xlabel,ylabel,xlim,ylim,grid,draw,plt
from matplotlib import animation
#from IPython.core.debugger import Tracer
#debug_here = Tracer()

L = 10
dims = 2
c = Container(dims,10)
initializer = ParticleInitialize()
c = initializer('line',c)

distance_matrix = PeroidicDistanceMatrix(L)
force = LeonardJonesForce(distance_matrix,c.masses)
forward = VerletIntegrator(.01)

num_frames = 1000
def next_frame(i):
  (dx,dv) = forward(force,c.x,c.v,i)
  c.integrate(dx,dv)
  c.draw(gca())
  return c.circles

circles = []
fig = plt.figure()
ax = plt.gca()
ax.set_aspect('equal')
anim = animation.FuncAnimation(fig,next_frame,frames=num_frames,blit=True)
plt.show()
#for i in range(5):
#  figure()
#  ax = gca()
#  ax.set_aspect('equal')
#  (dx,dv) = forward(force,c.x,c.v,i) 
#  c.integrate(dx,dv)
#  c.draw(ax)
#  show()

