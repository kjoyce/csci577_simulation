"""
Created on Sat Mar 02 15:14:08 2013

@author: Nathan Sponberg
"""
import csv
from numpy import array, hstack

class Container(object):
    def __init__(self, xdim, ydim, zdim = 0):
        self.Lx = xdim
        self.Ly = ydim
        self.Lz = zdim
        self.xpos = array([])
        self.ypos = array([])
        self.zpos = array([])
        self.xvel = array([])
        self.yvel = array([])
        self.zvel = array([])
        self.massVector = array([])
        self.numParticles = 0
        
    def addParticle(self,initState,mass):
        self.numParticles += 1
        self.xpos = hstack((self.xpos, initState[0]))
        self.ypos = hstack((self.ypos, initState[1]))
        self.zpos = hstack((self.zpos, initState[2]))
        self.xvel = hstack((self.xvel, initState[3]))
        self.yvel = hstack((self.yvel, initState[4]))
        self.zvel = hstack((self.zvel, initState[5]))
        self.massVector = hstack((self.massVector, mass))
    
    def summary(self):
        print self.xpos
        print self.ypos
        print self.zpos
        print self.xvel
        print self.yvel
        print self.zvel
        print self.massVector
        print self.numParticles
        print self.Lx
        print self.Ly
        print self.Lz
        
def InitalizeContainer(filename):
    with open(filename, 'r') as f:
      reader = csv.reader(f,delimiter=" ")
      numParticles = int(reader.next()[0])
      xdim = float(reader.next()[0])
      ydim = float(reader.next()[0])
      zdim = float(reader.next()[0])
      container = Container(xdim,ydim,zdim)
      for row in reader:
          state = [float(i) for i in row]
          container.addParticle(state[0:-1],state[-1])
      return container
        
test1 = InitalizeContainer("particle_test.txt")
test1.summary()  
