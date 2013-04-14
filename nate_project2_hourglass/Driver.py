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
num_frames = 15000
frameSkip = 2
neighborUpdateInterval = 1
data = array([0,0]) #junk variable, let in so as not to break code
count = 0
rad = (2.**(1./6.))/2.
Lx = 20.
Ly = 20.

floorSize = 18
wallSize = 3
slantSize = 10

###Use Radians!!!
slantStart = 20.
slantDegree = pi/3.1
slantMultiplierX = 2*rad*cos(slantDegree)
slantMultiplierY = 2*rad*sin(slantDegree)

container = Container(floorSize,wallSize,slantSize,Lx,Ly)
dist = container.Lx / 5.
vel = dist /5.
animate = True
iterationTimed = 500
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
container.openingPosition =2




for i in range(floorSize - 7):
    container.addParticle((i*2*(rad) + 4*rad, slantStart,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
for i in range(floorSize - 8):
    container.addParticle((i*2*(rad) + 5*rad, slantStart - 2*rad,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
for i in range(floorSize - 9):
    container.addParticle((i*2*(rad) + 6*rad, slantStart - 4*rad,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)

for i in range(floorSize - 10):
    container.addParticle((i*2*(rad) + 7*rad, slantStart - 6*rad,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)
for i in range(floorSize - 11):
    container.addParticle((i*2*(rad) + 8*rad, slantStart - 8*rad,0.,
                          0.,0.,0.,
                          0.,0.,0.),
                          1.)

force = GranularForces()

#xacl,yacl,zacl = force(container)
#container.xacl = xacl
#container.yacl = yacl
#container.zacl = zacl

integrator = Verlet(dt)

setNeighbors = Neighbors()
setNeighbors.UpdateNeighbors(container, cutOff)
    

if animate == True:

    circles = []
    fig = plt.figure(figsize = (10,10))
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
    
  
    
    def next_frame(i, data, count): 
        ##integration loop
      #for i in range(frameSkip):
          #integrator(force,container)
          #print container.neighborList
      for j in range(frameSkip):
          count += 1
          if count%neighborUpdateInterval == 0:
              setNeighbors.UpdateNeighbors(container, cutOff)
          integrator(force,container)
      for i in range(len(circles)):
        circles[i].center = (container.xpos[i], container.ypos[i])
      return circles
    
      #debug_here()
    anim = animation.FuncAnimation(fig,next_frame,init_func=init,
                                   frames=num_frames, fargs=(data,count), interval=1
                                   ,blit=True)
    
    plt.show()

else:
    start = time()
    setNeighbors = Neighbors()
    setNeighbors.UpdateNeighbors(container)
    for i in range(iterationTimed):    
        if i%neighborUpdateInterval == 0:
          setNeighbors.UpdateNeighbors(container)
        integrator(force,container)
        
    elapsed = time() - start
    print elapsed
