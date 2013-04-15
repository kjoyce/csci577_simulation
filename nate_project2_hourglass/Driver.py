# -*- coding: utf-8 -*-
"""
Created on Tue Mar 05 19:27:31 2013

@author: Nathan Sponberg
"""

from Integrate import Verlet
from Force import GranularForces
from Container import InitalizeContainer, Container
from DistMatrix import Neighbors
import matplotlib as plt
from matplotlib import animation
from pylab import * # for plotting commands
from time import time
from scipy import sin,cos,pi
from IPython.core.debugger import Tracer
debug_here = Tracer()

##########################
####Initial Conditions####
END_TIME = 1.
t = 0.
dt = 0.010
initConditions = "square_lattice"
num_frames = 100
frameSkip = 1
neighborUpdateInterval = 1
data = array([0,0]) #junk variable, let in so as not to break code
count = 0
rad = (2.**(1./6.))/2.
Lx = 27./sqrt(3)
Ly = 20.

floorSize = 18
wallSize = 3
slantSize = 10

###Use Radians!!!
slantStart = 20.
slantDegree = pi/3.1
slantMultiplierX = 2*rad*cos(slantDegree)
slantMultiplierY = 2*rad*sin(slantDegree)
gamma = 20

num_rows = 6

container = Container(floorSize,wallSize,slantSize,Lx,Ly)
dist = container.Lx / 5.
vel = dist /5.
animate = False
iterationTimed = 1450
cutOff = 1.*(2*rad)
#########################
"""
for i in range(floorSize):
    container.addParticle((i*2*(rad) + rad/2, 1.,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)

for i in range(wallSize):
    container.addParticle((rad/2, 1.+(2*rad*(i+1)),0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
for i in range(wallSize):
    container.addParticle((rad/2 + (floorSize-1)*2*rad, 1.+(2*rad*(i+1)),0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
"""
'''                        
for i in range(slantSize):
    container.addParticle((rad/2 +(i*slantMultiplierX), slantStart-(i*slantMultiplierY),0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
for i in range(slantSize):
    container.addParticle((rad/2 -(i*slantMultiplierX)+(floorSize-6)*2*rad, slantStart-(i*slantMultiplierY),0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
    #if i == (slantSize - 1):
        ####this doesn't work right now
        #container.openingWidth = container.xpos[-1] - container.xpos[-2]
        ################
'''
container.openingPosition = 8




def set_particles(rows):
  for j in range(rows):
    for i in range(floorSize - (7+j)):
	container.addParticle((i*2*(rad) + (4+j)*rad, slantStart - 2*j*rad,0.,
			      0.,0.,0.,
			      0.,0.,0.),
			      1.)

set_particles(num_rows)
force = GranularForces(gamma=gamma)

#xacl,yacl,zacl = force(container)
#container.xacl = xacl
#container.yacl = yacl
#container.zacl = zacl

integrator = Verlet(dt)

setNeighbors = Neighbors()
setNeighbors.UpdateNeighbors(container, cutOff)
    

if animate == True:

    circles = []
    fig = plt.figure(figsize = (10,8))
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim((0,container.Lx))
    ax.set_ylim((0,container.Ly))
    ax.plot([0,10/sqrt(3)],[20,10],'b-')
    ax.plot([17/sqrt(3),27/sqrt(3)],[10,20],'b-')
    
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
    
    for i in range(container.numParticles):
      e = Circle( (container.xpos[i],container.ypos[i]), radius=rad)#, animated = True)
      e = prettify_circle(e)
      circles.append(ax.add_patch(e))
    def init():
      return circles
    
    def next_frame(i): 
        ##integration loop
      #for i in range(frameSkip):
          #integrator(force,container)
          #print container.neighborList
      for j in range(frameSkip):
	  global count
          count += 1
	  if i % 100 == 0:
	    print "frame: {}".format(count)
          if count%neighborUpdateInterval == 0:
              setNeighbors.UpdateNeighbors(container, cutOff)
          integrator(force,container)
	  if count % 200 == 0:
	    title("Frame: {}".format(count))
	    savefig("hourglass_{}.pdf".format(count))
      for i in range(len(circles)):
        circles[i].center = (container.xpos[i], container.ypos[i])
      return circles
    
      #debug_here()
    anim = animation.FuncAnimation(fig,next_frame,init_func=init,
                                   frames=num_frames,interval=1)
#                                   ,blit=True)
    
    plt.show()

else:
    bridge_height_estimate = zeros((iterationTimed,10))
    for i in range(iterationTimed):    
        #if i%neighborUpdateInterval == 0:
	setNeighbors.UpdateNeighbors(container,cutOff)
        integrator(force,container)
	xslice_acc = []
	yslice_acc = []
	for h in range(10,20):
	  idx = (container.ypos < h+1)*(container.ypos >= h)
	  yslice_acc.append(average(container.yacl[idx]))
	  xslice_acc.append(average(container.xacl[idx]))
	xslice_acc = array(xslice_acc)
	yslice_acc = array(yslice_acc)
	xslice_acc[isnan(xslice_acc)] = 0
	yslice_acc[isnan(yslice_acc)] = 0
	bridge_height_estimate[i] = yslice_acc 
	if i % 100 == 0 and i > 100:
	  print "frame: {}".format(i)
	  f = figure(figsize=(9,5))
	  subplot(121)
	  plot(range(10,20),xslice_acc)
	  title("Average Hor. Acceleration vs. Height ")
	  xlabel("Height")
	  ylabel("Average X-Acceleration")
	  ylim(-5,15)
	  subplot(122)
	  plot(range(10,20),yslice_acc)
	  title("Average Ver. Acceleration vs. Height ")
	  xlabel("Height")
	  ylabel("Average Y-Acceleration")
	  subplots_adjust(wspace=.4)  # Note this makes space
	  ylim(-5,15)
	  savefig("component_slice_forces_{}.pdf".format(i))

    figure()
    plot(container.particleFlux)
    ylim(-1,5)
    title("Granular Flow Rate") 
    xlabel("Time")
    ylabel("Particle Flux through $y=8$ to $y=10$") 
    savefig("particle_flux.pdf")

    import pickle
    f = open('yaccel.dump','w')
    pickle.dump([bridge_height_estimate,container.particleFlux],f)
    f.close()
# After unpickeling
#In [105]: contourf(arange(1450),arange(10,20),tot_yaccel.T) 
#Out[105]: <matplotlib.contour.QuadContourSet instance at 0xbb2d20c>
#
#In [106]: colorbar() 
#Out[106]: <matplotlib.colorbar.Colorbar instance at 0xbc7058c>
#
#In [107]: title("Avg. Vertical of Force as Various Heights") 
#Out[107]: <matplotlib.text.Text at 0xb741a8c>
#
#In [108]: ylabel("Height") 
#Out[108]: <matplotlib.text.Text at 0xbd2d3cc>
#
#In [109]: xlabel("Time") 
#Out[109]: <matplotlib.text.Text at 0xba5388c>
#
#In [110]: savefig("jam_show.pdf") 

