from Initialize import ParticleInitialize
from pylab import figure,clf,plot,title,show,gca,xlabel,ylabel,xlim,ylim,grid,draw,plt,Circle,linspace,zeros,arange,ones,ogrid,nditer,floor,ceil,Rectangle
from matplotlib import animation
#from IPython.core.debugger import Tracer
#debug_here = Tracer()

########## PARAMS ###############
num_frames = 1000
ekg_length = 80
delay      = 0
vid_format = ".mp4"

run_backwards  = False
save_animation = True
print_frame    = True
tile_domain    = True
crunch_domain  = False

#case = "one"
#case = "two"
#case = "three"
case = "four"
#case = "six"
#case = "eight"
#case = "line"
#case = "square_lattice"
#case = "triangle_lattice"
#################################
initializer = ParticleInitialize()
c,distance_matrix,force,integrate,xlim,ylim,pot_energy_lim,kin_energy_lim,tot_energy_lim,pressure_lim = initializer(case)
print pot_energy_lim
print kin_energy_lim
print tot_energy_lim
Linit = c.L.copy()
file_name = case

########### ANIMATION  ##########
circles = []
fig = plt.figure(figsize=(13,6))
ax = plt.subplot2grid((2,4),(0,0),colspan=2,rowspan=2,xlim=xlim,ylim=ylim,aspect='equal')
a_title = ax.set_title("",animated=not(save_animation))
ax2 = plt.subplot2grid((2,4),(0,2),title='Potential Energy',ylim=pot_energy_lim)
ax3 = plt.subplot2grid((2,4),(1,2),title='Kinetic Energy',ylim=kin_energy_lim)
ax4 = plt.subplot2grid((2,4),(0,3),title='Total Energy',ylim=tot_energy_lim)
ax5 = plt.subplot2grid((2,4),(1,3),title='Pressure',ylim=pot_energy_lim)

for a in ax2,ax3,ax4,ax5:
  a.grid()
  a.set_xticks([])
fig.subplots_adjust(wspace=.4)  # this makes space between subplots
fig.subplots_adjust(hspace=.4)  # this makes space between subplots

# Plot a square around the primary calculation domain
domain = ax.add_patch(Rectangle((0,0),c.L[0],c.L[1], fc="none", alpha=.6, ec="lightsteelblue",animated=True))

def my_circle(x,y):
  e = Circle( (x,y), radius=2**(-.5/.6), animated=not(save_animation), clip_on=True)
  color="lightsteelblue"
  facecolor="green"
  alpha=.6
  e.set_clip_box(ax.bbox)
  e.set_edgecolor( color )
  e.set_linewidth(3)
  e.set_facecolor( facecolor )  # "none" not None
  e.set_alpha( alpha )
  return e

def circle_iter_repeat(vec,xlim,ylim,LL):
  """This will iterate through all possible circles within xlim + (-1,1) and ylim + (-1,1)"""
  left   = (int(xlim[0]/LL[0])-1)*LL[0]  
  right  = (int(xlim[1]/LL[0])+1)*LL[0]
  bottom = (int(ylim[0]/LL[1])-1)*LL[1]  
  top    = (int(ylim[1]/LL[1])+1)*LL[1]   
  coords = ogrid[left:right:LL[0],bottom:top:LL[1]]
  i = 0      
  for x,y in nditer(coords):
    for dx,dy in (vec[:,0:2]): 
      yield x+dx,y+dy,i
      i += 1

def circle_iter_single(vec,xlim,ylim,LL):
  for i in range(len(vec)):
    yield vec[i,0],vec[i,1],i

if tile_domain:
  circle_iter = circle_iter_repeat
else:
  circle_iter = circle_iter_single

for x,y,i in circle_iter(c.x,xlim,ylim,c.L):
  e = my_circle(x,y)
  circles.append(ax.add_patch(e))
  
print "Num Circles: {}".format(len(circles))
# Set up energy lines
j = 0
t = arange(ekg_length)
potential_dat    = zeros(t.shape)
kinetic_dat      = zeros(t.shape) 
total_energy_dat = zeros(t.shape) 
pressure_dat     = zeros(t.shape) 
lines = []
lines.append(ax2.plot(t,potential_dat,'b-',animated=not(save_animation))[0])
lines.append(ax3.plot(t,kinetic_dat,'b-',animated=not(save_animation))[0])
lines.append(ax4.plot(t,total_energy_dat,'b-',animated=not(save_animation))[0])
lines.append(ax5.plot(t,pressure_dat,'b-',animated=not(save_animation))[0])

def next_frame(i):
  global j
  backwards = run_backwards and i > num_frames/2
  start_backwards = i == int(num_frames/2) and run_backwards
  if backwards:
    direction = "backward"
    (dx,dv) = integrate.backward(force,c.x,c.v,i)
  else: 
    direction = "forward"
    (dx,dv) = integrate.forward(force,c.x,c.v,i)
  if start_backwards:
    for e in circles:
      e.set_facecolor("red")
    for l in lines:
      l.set_color("red")
  if i == 0:
    for e in circles:
      e.set_facecolor("green")
    for l in lines:
      l.set_color("green")
  c.integrate(dx,dv)
  if print_frame:
    print "Frame: {} {}".format(i,direction)
  line_dir = (-1)**backwards
  potential_dat[(line_dir*j)%ekg_length] = force.potential_energy
  kinetic_dat  [(line_dir*j)%ekg_length] = force.kinetic_energy
  pressure_dat [(line_dir*j)%ekg_length] = force.pressure
  total_energy_dat = potential_dat + kinetic_dat
  j = j + 1
  lines[0].set_data(t,potential_dat    [(t+line_dir*j)%ekg_length])
  lines[1].set_data(t,kinetic_dat      [(t+line_dir*j)%ekg_length])
  lines[2].set_data(t,total_energy_dat [(t+line_dir*j)%ekg_length])
  lines[3].set_data(t,pressure_dat     [(t+line_dir*j)%ekg_length])
  a_title.set_text("Frame: {}".format(j))
  if crunch_domain:
    crunch_rate = .999
    c.L *= crunch_rate
    force.distance_matrix.L *= crunch_rate
    domain.set_width (c.L[0])
    domain.set_height(c.L[1])
  for x,y,i in circle_iter(c.x,xlim,ylim,Linit):
    circles[i].center = (x,y)
  return tuple(circles + lines +[a_title,domain])
 
#next_frame(0)
anim = animation.FuncAnimation(fig,next_frame,frames=num_frames,interval=delay,blit=True)

if save_animation:
  anim.save(file_name+vid_format,fps=45)
else:
  plt.show()

