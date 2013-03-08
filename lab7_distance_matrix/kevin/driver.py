from Container import Container
from Force import LeonardJonesForce
from Integrator import VerletIntegrator
from DistanceMatrix import PeroidicDistanceMatrix
from Initialize import ParticleInitialize
from pylab import figure,clf,plot,title,show,gca,xlabel,ylabel,xlim,ylim,grid,draw,plt,Circle,linspace,zeros,arange,ones
from matplotlib import animation
#from IPython.core.debugger import Tracer
#debug_here = Tracer()

########## PARAMS ###############
num_frames = 800
potential_ylim = 250
ekg_length = 25
run_backwards = False
delay = 1
save_animation = False
print_frame = False
file_name = "square_lattice"
L = 10
dims = 2
initializer = ParticleInitialize()
c = Container(dims,10)
c = initializer(file_name,c)
#################################

distance_matrix = PeroidicDistanceMatrix(L)
force = LeonardJonesForce(distance_matrix,c.masses)
integrate = VerletIntegrator(.01)
########### FORCED ANIMATION ##########
circles = []
fig = plt.figure(figsize=(10,5))
ax = plt.subplot2grid((2,2),(0,0),rowspan=2,xlim=(0,c.L[0]),ylim=(0,c.L[1]),aspect='equal')
ax2 = plt.subplot2grid((2,2),(0,1),title='Potential Energy',ylim=(-1,potential_ylim))
ax2.grid()
ax2.set_xticklabels([])
ax3 = plt.subplot2grid((2,2),(1,1),title='Kinetic Energy')
ax3.grid()
ax3.set_xticklabels([])
ax3.text(.4,.5,"????",fontsize=25, transform=ax3.transAxes)
fig.subplots_adjust(wspace=.4)  # this makes space between subplots
fig.subplots_adjust(hspace=.4)  # this makes space between subplots

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

# Pre initializing is necessary I guess
for x in c.x:
  e = Circle( (x[0],x[1]), radius=2**(-.5/.6), animated=not(save_animation))
  e = prettify_circle(e)
  circles.append(ax.add_patch(e))

# set up energy line
j = 0
t = arange(ekg_length)
potential_dat = zeros(t.shape)
#line = (ax2.plot([],[],'b-', animated=not(save_animation))[0])
line = (ax2.plot(t,potential_dat*300,'b.',animated=not(save_animation))[0])

def init():
  return tuple(circles + [line])


def next_frame(i):
  global j
  backwards = run_backwards and i>= num_frames/2
  if backwards:
    direction = "backward"
    (dx,dv) = integrate.backward(force,c.x,c.v,i)
  else: 
    direction = "forward"
    (dx,dv) = integrate.forward(force,c.x,c.v,i)
  c.integrate(dx,dv)
  if print_frame:
    print "Frame: {} {}".format(i,direction)
  for i in range(len(circles)):
    circles[i].center = (c.x[i,0],c.x[i,1])
    if backwards:
      circles[i].set_facecolor("red")
    else:
      circles[i].set_facecolor("green")
  # update potential
  potential_dat[(j%ekg_length)] = force.potential_energy
  j = j + 1
  line.set_data((t+j)%ekg_length,potential_dat)
  return tuple(circles + [line])
 
anim = animation.FuncAnimation(fig,next_frame,init_func=init,frames=num_frames,interval=delay,blit=True)

if save_animation:
  anim.save(file_name+".mp4",fps=45)
else:
  plt.show()

