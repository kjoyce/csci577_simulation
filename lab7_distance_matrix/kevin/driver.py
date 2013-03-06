from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from Initialize import ParticleInitialize
from pylab import figure,clf,plot,title,show,gca,xlabel,ylabel,xlim,ylim,grid,draw,plt,Circle
from matplotlib import animation
from IPython.core.debugger import Tracer
debug_here = Tracer()

L = 10
dims = 2
c = Container(dims,10)
initializer = ParticleInitialize()
c = initializer('line',c)

distance_matrix = PeroidicDistanceMatrix(L)
force = LeonardJonesForce(distance_matrix,c.masses)
forward = VerletIntegrator(.01)

########### FORCED ANIMATION ##########
########## PARAMS ###############
num_frames = 350
delay = 1
save_animation = False
file_name = "line"
#################################

circles = []
fig = plt.figure()
ax = plt.gca()
ax.set_aspect('equal')
ax.set_xlim((0,c.L[0]))
ax.set_ylim((0,c.L[1]))

def prettify_circle(e):
  color="lightsteelblue"
  facecolor="green"
  alpha=.6
  e.set_clip_box(ax.bbox)
  e.set_edgecolor( color )
  e.set_linewidth(3)
  e.set_facecolor( facecolor )  # "none" not None
  e.set_alpha( alpha )
  return e

## Pre initializing is necessary I guess
for x in c.x:
  e = Circle( (x[1],x[0]), radius=2**.2, animated=True)
  e = prettify_circle(e)
  circles.append(ax.add_patch(e))
def init():
  return circles

def next_frame(i):
  print "Frame: {}".format(i)
  (dx,dv) = forward(force,c.x,c.v,i)
  c.integrate(dx,dv)
  for i in range(len(circles)):
    circles[i].center = (c.x[i,1],c.x[i,0])
  return circles
  
anim = animation.FuncAnimation(fig,next_frame,init_func=init,frames=num_frames,interval=delay,blit=True)

if save_animation:
  anim.save(file_name+".mp4",fps=25)
else:
  plt.show()

