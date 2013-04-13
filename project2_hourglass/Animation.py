from numpy import linspace,zeros,arange,ones,ogrid,nditer,floor,ceil
import matplotlib.pyplot as plt
import matplotlib.animation as animation
class Container_Animation(object):
  def __init__(self,container,integrate,force,xlim,ylim,
	       pull_force_lim,ave_vel_lim,
	       num_frames=1000,ekg_length=80,delay=0,
	       save_animation=False,run_backwards=False,
	       print_frame=False,tile_domain=False,
	       crunch_domain=False,vid_format="mp4",filename="out"):
    self.fig = plt.figure(figsize=(13,6))
    self.c = container
    self.integrate = integrate
    self.force = force
    self.circles = []
    self.axes = self._set_axes(xlim,ylim,pull_force_lim,ave_vel_lim)
    self.ekg_length = ekg_length
    self.save_animation = save_animation
    self.lines = self._set_lines()
    self.run_backwards = run_backwards
    self.num_frames = num_frames
    self.print_frame = print_frame
    self.j = 0
    self.crunch_domain = crunch_domain
    self.xlim = xlim
    self.ylim = ylim
    self.Linit = self.c.L.copy()
    self.title = self.axes[0].set_title("",animated=not(save_animation))
    self.tile_domain = tile_domain
    self.delay = delay
    # Plot a square around the primary calculation domain
    self.domain = self.axes[0].add_patch(plt.Rectangle((0,0),self.c.L[0],self.c.L[1], fc="none", alpha=.6, ec="lightsteelblue",animated=not(save_animation)))
    if tile_domain:
      self.circle_iter = self.circle_iter_repeat
    else:
      self.circle_iter = self.circle_iter_single
    for x,y,i in self.circle_iter(self.c.x,xlim,ylim,self.c.L):
      e = self.make_circle(x,y)
      self.circles.append(self.axes[0].add_patch(e))
    self.j = 0
    self.filename = filename
    self.vid_format = vid_format
    self.anim = animation.FuncAnimation(self.fig,self._animate_frame,frames=self.num_frames,blit=True)

  def _set_axes(self,xlim,ylim,pull_force_lim,ave_vel_lim):
    axes = []
    axes.append(plt.subplot2grid((2,3),(0,0),colspan=2,rowspan=2,xlim=xlim,ylim=ylim,aspect='equal'))
    axes.append(plt.subplot2grid((2,3),(0,2),title='Pulling Force',ylim=pull_force_lim))
    axes.append(plt.subplot2grid((2,3),(1,2),title='Average Velocity',ylim=ave_vel_lim))
    for a in axes[1:]:
      a.grid()
      a.set_xticks([])
    return axes

  def _set_lines(self):
    ekg_length = self.ekg_length
    t = arange(ekg_length)
    self.t = t
    self.pull_force_dat    = zeros(t.shape)
    self.ave_vel_dat       = zeros(t.shape) 
    lines = []
    axes = self.axes
    save_animation = self.save_animation
    lines.append(axes[1].plot(t,self.pull_force_dat,'b-',animated=not(save_animation))[0])
    lines.append(axes[2].plot(t,self.ave_vel_dat,'b-',animated=not(save_animation))[0])
#    lines.append(axes[0].plot([],[],'-',lw=2,animated=not(save_animation))[0])
    return lines
    
  def circle_iter_repeat(self,vec,xlim,ylim,LL):
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

  def circle_iter_single(self,vec,xlim,ylim,LL):
    for i in range(len(vec)):
      yield vec[i,0],vec[i,1],i
    
  def make_circle(self,x,y):
    e = plt.Circle( (x,y), radius=2**(-.5/.6), animated=not(self.save_animation), clip_on=True)
    color="lightsteelblue"
    facecolor="green"
    alpha=.6
    e.set_clip_box(self.axes[0].bbox)
    e.set_edgecolor( color )
    e.set_linewidth(3)
    e.set_facecolor( facecolor )  # "none" not None
    e.set_alpha( alpha )
    return e

  def _update_lines(self,backwards):
    lines,ekg_length,force,j,pull_force_dat,ave_vel_dat,t = self.lines,self.ekg_length,self.force,self.j,self.pull_force_dat,self.ave_vel_dat,self.t
    line_dir = (-1)**backwards
    pull_force_dat[(line_dir*j)%ekg_length] = force.pull_force[-1]
    ave_vel_dat[(line_dir*j)%ekg_length] = self.c.avg_sled_velocity
    j = j + 1
    lines[0].set_data(t,pull_force_dat[(t+line_dir*j)%ekg_length])
    lines[1].set_data(t,ave_vel_dat   [(t+line_dir*j)%ekg_length])
#    lines[2].set_data([self.c.x[-1,0],self.c.t*.1+self.force.xinit],[self.c.x[-1,1],1])

  def _animate_frame(self,i):
    self.j += 1
    run_backwards,num_frames,integrate,force,domain,circles,c,j,crunch_domain,xlim,ylim,Linit,lines = self.run_backwards,self.num_frames,self.integrate,self.force,self.domain,self.circles,self.c,self.j,self.crunch_domain,self.xlim,self.ylim,self.Linit,self.lines
    backwards = run_backwards and i > num_frames/2
    start_backwards = i == int(num_frames/2) and run_backwards
    if backwards:
      direction = "backward"
      c.etargetni()
    else: 
      direction = "forward"
      c.integrate()
    if start_backwards:
      for e in circles:
	e.set_facecolor("red")
      for l in lines:
	l.set_color("red")
    if i == 0:
      for e in circles:
	e.set_facecolor("green")
      for l in self.lines:
	l.set_color("green")
    self._update_lines(backwards)
    self.title.set_text("Frame: {}".format(j))
    if crunch_domain and not(j%20):
      crunch_rate = .995
      c.updateL(c.L*crunch_rate)
      domain.set_width (c.L[0])
      domain.set_height(c.L[1])
    if self.tile_domain:
      circle_iter = self.circle_iter_repeat
    else:
      circle_iter = self.circle_iter_single
    for x,y,i in circle_iter(c.x,xlim,ylim,Linit):
      circles[i].center = (x,y)
      if i == c.hot_idx:  # this is a dirty hack that doesn't work when the domain is tiled
	circles[i].set_facecolor("red")
    if self.print_frame:
      print "Frame: {} {} {}".format(j,direction,num_frames)
      print "PE: {} KE: {} TE: {} P: {}".format(force.potential_energy,force.kinetic_energy,force.potential_energy+force.kinetic_energy,force.pressure)
    return tuple(circles + self.lines +[self.title,domain])

  def show(self):
    plt.show()

  def save(self):
    self.anim.save(self.filename+self.vid_format,fps=45)
