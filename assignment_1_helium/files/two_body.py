import csv
from scipy.integrate import odeint  # for integrate.odeint
from pylab import *  # for plotting commands
from IPython.core.debugger import Tracer
from matplotlib.pyplot import *
from matplotlib import animation
debug_here = Tracer()

############# Integrating parameters ################
num_states = 5000
time_length = 89
#file_list = ['helium2.txt','helium3.txt','helium4.txt','helium1.txt']
#file_list = ['helium3a.txt']#,'helium2.txt','helium3.txt','helium4.txt']
#file_list = ['scrambled_egg.txt']
#file_list = ['third_times_a_charm.txt']
file_list = ['butterfly.txt']

####### Animation Parameters ###########
save_animation = True
figdim = (13,6)  # dimensions of figure
frameskip = 10	 # animate every nth frame
delay = 1        # delay between redrawing frames in miliseconds
percent_cap = 1  # This is what percentage of data to capture
margin = .1	 # this is how much of a margin from outlier (in terms of spread)
max_xlim = -7,7
max_ylim = -7,7

#######################################
def two_body(state,t):  #OUCH! The signature is reversed for odeint!
  x1 = array([state[0],state[1]])  # body 1 position
  v1 = array([state[2],state[3]])  # body 1 velocity
  x2 = array([state[4],state[5]])
  v2 = array([state[6],state[7]])

  r1  = norm(x1)
  r2  = norm(x2)
  r21 = norm(x2 - x1)

  a1 = -2.*x1/r1**3 + (x1-x2)/r21**3
  a2 = -2.*x2/r2**3 + (x2-x1)/r21**3
  derivatives = array([v1,a1,v2,a2]).T.flatten(1)
  return derivatives
 
def read_file(filename):
  with open(filename, 'r') as f:
    reader = csv.reader(f,delimiter=' ')
    numBodies = int(reader.next()[0]) # maybe do some check later but useless as is
    bod = []
    for row in reader:
      temp = [float(i) for i in row]  
      bod.append(temp[0:-1])
    bod = array(bod)
    return bod
for fname in file_list:
  xinit = read_file(fname)
  times = linspace(0,time_length,num_states)
  xx = odeint(two_body,xinit.flatten(),times).T #,rtol=10**-13,atol=10**-13).T

  r1 = sqrt(xx[0]**2 + xx[1]**2)
  r2 = sqrt(xx[4]**2 + xx[5]**2)
  r21 = sqrt( (xx[4]-xx[0])**2 + (xx[5]-xx[1])**2) 
  E = .5*((xx[2]**2 + xx[3]**2) + (xx[6]**2+xx[7]**2)) + ( -2/r1  -2/r2 + 1/r21 ) 
  if E[0] == 0:
    deltaE = E
    debug_here()
  else:
    deltaE = (E-E[0])/E[0]

  L = (xx[0]*xx[3] - xx[1]*xx[2]) + (xx[4]*xx[7] - xx[5]*xx[6])
  if L[0] == 0:
    deltaL = L
  else:
    deltaL = (L-L[0])/L[0]

################ Animation ##################
  points = []
  lines = []
  fig = plt.figure(figsize = figdim)
  ax1 = plt.subplot2grid((2,2),(0,1))
  ax1.plot(times,deltaE,'b:')
  xlim(times[0],times[-1])
  xlabel('Time')
  ylabel('% $\Delta E /M$')
  grid()

  ax2 = plt.subplot2grid((2,2),(1,1))
  ax2.plot(times,deltaL,'b:')
  xlabel('Time')
  ylabel('% $\Delta L / M$')
  grid()

  subplots_adjust(wspace=.4)  # this makes space between subplots
  subplots_adjust(hspace=.4)  # this makes space between subplots

  idx = floor((1-percent_cap)*size(times))
  xbnds = sort(hstack([xx[0],xx[4]]))
  sp = xbnds[-idx-1]-xbnds[idx]
  xlims = -margin*sp+xbnds[idx],margin*sp+xbnds[-idx-1]
  ybnds = sort(hstack([xx[1],xx[5]]))
  sp = ybnds[-idx-1]-ybnds[idx]
  ylims = -margin*sp+ybnds[idx],margin*sp+ybnds[-idx-1]
  if ( xlims[0] < max_xlim[0]\
    or ylims[0] < max_ylim[0]\
    or xlims[1] > max_xlim[1]\
    or ylims[1] > max_ylim[1]):
    xlims = max_xlim
    ylims = max_ylim
  ax = plt.subplot2grid((2,2),(0,0),rowspan=2,xlim=xlims, ylim=ylims)
  points.append(ax.plot([], [], 'r.')[0]) # initialize red electron
  points.append(ax.plot([], [], 'b.')[0]) # initialize blue electron

  lines.append(ax.plot([], [],'r-')[0]) # initialize red electron trail
  lines.append(ax.plot([], [],'b-')[0]) # initialize blue electron trail
  lines.append(ax1.plot([], [],'b-')[0]) # initialize energy line
  lines.append(ax2.plot([], [],'b-')[0]) # initialize momentum line

  ax.plot(0,0,'go',markersize=20) # Plot nucleus
  title('Helium Atom Model')
  xlabel('Horizontal Distance')
  ylabel('Vertical Distance')


# Initial animation state
  def init():
    for line in lines:
      line.set_data([],[])
    for point in points:
      point.set_data([], [])
    return tuple(lines + points)

  def animateTrace(i):
      idx = frameskip*i
      lines[0].set_data(xx[0,0:idx], xx[1,0:idx]) # update e1 trail
      lines[1].set_data(xx[4,0:idx], xx[5,0:idx]) # update e2 trail
      lines[2].set_data(times[0:idx],deltaE[0:idx])  # update energy trail
      lines[3].set_data(times[0:idx],deltaL[0:idx])  # update momentum trail
	  
      points[0].set_data(xx[0,idx],xx[1,idx]) # update e1 pt
      points[1].set_data(xx[4,idx],xx[5,idx]) # update e2 pt
      return tuple(lines + points)

  anim = animation.FuncAnimation(fig, animateTrace, init_func=init, frames=size(times)/frameskip, interval=delay, blit=True) 
  fname = fname.rstrip('.txt')  
  print fname
  if save_animation:
    anim.save(fname+'.mp4',fps=25)#,extra_args=['-vcodec', 'libx264'])
  else: 
    plt.show()
