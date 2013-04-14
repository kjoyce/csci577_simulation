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
gamma = 30

num_rows = 6

container = Container(floorSize,wallSize,slantSize,Lx,Ly)
dist = container.Lx / 5.
vel = dist /5.
animate = False
iterationTimed = 1200
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
      e = Circle( (container.xpos[i],container.ypos[i]), radius=rad, animated = True)
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
      for i in range(len(circles)):
        circles[i].center = (container.xpos[i], container.ypos[i])
      return circles
    
      #debug_here()
    anim = animation.FuncAnimation(fig,next_frame,init_func=init,
                                   frames=num_frames,interval=1
                                   ,blit=True)
    
    plt.show()

else:
    bridge_height_estimate = array([])
    for i in range(iterationTimed):    
        #if i%neighborUpdateInterval == 0:
	setNeighbors.UpdateNeighbors(container,cutOff)
        integrator(force,container)
	print "frame: {}".format(i)
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
	idx = find(yslice_acc > 0)
	if not(idx is None):
	  bridge_height_estimate = hstack((bridge_height_estimate, 10+idx) )
	if i % 100 == 0 and i > 100:
	  figure()
	  subplot(121)
	  plot(range(10,20),xslice_acc)
	  title("Average Horizontal Acceleration vs. Height ")
	  xlabel("Height")
	  ylabel("Average X-Acceleration")
	  subplot(122)
	  plot(range(10,20),yslice_acc)
	  title("Average Vertical Acceleration vs. Height ")
	  xlabel("Height")
	  ylabel("Average Y-Acceleration")
	  savefig("component_slice_forces_{}.pdf".format(i))

    figure()
    plot(container.particleFlux)

    figure()
    hist(bridge_height_estimate,bins=8)

        
