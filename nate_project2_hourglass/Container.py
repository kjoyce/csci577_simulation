# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 15:14:08 2013

@author: Nathan Sponberg
"""
import csv
from numpy import array, hstack, size, tile, ones

class Container(object):
    def __init__(self, floorSize, wallSize, slantSize, xdim, ydim, zdim = 0, springEqDist = 2.*2.**(1./6.)):
        self.floorSize = floorSize
        self.wallSize = wallSize
        self.slantSize = slantSize
        self.Lx = xdim
        self.Ly = ydim
        self.Lz = zdim
        self.springEqDist = springEqDist #equilibrium distance between lattice particles
        self.initDragPos = 0. #initial position of right most lattice particle
        ###postions
        self.xpos = array([],dtype='float64')
        self.ypos = array([],dtype='float64')
        self.zpos = array([],dtype='float64')
        ###velocities
        self.xvel = array([],dtype='float64')
        self.yvel = array([],dtype='float64')
        self.zvel = array([],dtype='float64')
        ###accelerations
        self.xacl = array([],dtype='float64')
        self.yacl = array([],dtype='float64')
        self.zacl = array([],dtype='float64')
        ###initial system conditions, used to reset container for multiple simulations
        self.xposinit = array([],dtype='float64')
        self.yposinit = array([],dtype='float64')
        self.zposinit = array([],dtype='float64')
        self.xvelinit = array([],dtype='float64')
        self.yvelinit = array([],dtype='float64')
        self.zvelinit = array([],dtype='float64')
        self.xaclinit = array([],dtype='float64')
        self.yaclinit = array([],dtype='float64')
        self.zaclinit = array([],dtype='float64')
        ###masses
        self.massVector = array([],dtype='float64')
        self.numParticles = 0
        ###empty variable used later
        self.adjMatrix = array([],dtype='bool')
        self.dragForce = []
        self.maxDragForces = []
        ###
        self.openingWidth = 0.
        self.openingPosition = 0.
        self.particleFlux = list()
        self.integrationIteration = 0

    def addParticle(self,initState,mass): #add particles to the system
        self.numParticles += 1
        self.xpos = hstack((self.xpos, initState[0]))
        self.ypos = hstack((self.ypos, initState[1]))
        self.zpos = hstack((self.zpos, initState[2]))
        self.xvel = hstack((self.xvel, initState[3]))
        self.yvel = hstack((self.yvel, initState[4]))
        self.zvel = hstack((self.zvel, initState[5]))
        self.xacl = hstack((self.xacl, initState[6]))
        self.yacl = hstack((self.yacl, initState[7]))
        self.zacl = hstack((self.zacl, initState[8]))
        self.xposinit = self.xpos
        self.yposinit = self.ypos
        self.zposinit = self.zpos
        self.xvelinit = self.xvel
        self.yvelinit = self.yvel
        self.zvelinit = self.zvel
        self.xaclinit = self.xacl
        self.yaclinit = self.yacl
        self.zaclinit = self.zacl

        self.massVector = hstack((self.massVector, mass))

    def getMassMatrix(self): ###not used any more, possibly useful in the future
        massMatrix = tile(self.massVector, (size(self.massVector),1))
        massMatrix = massMatrix.T
        return massMatrix

    ###sets adjacency martix for particle lattice, used to compute which particles
    ###have stiff springs connecting them
    def setAdjacency(self): 
        self.initDragPos = self.xpos[-1]
        self.adjMatrix = ones((self.numParticles,self.numParticles), dtype = bool)
        ###the order in which lattice particles are added (see driver class) means that for a given
        ###paticle in the lattice is connected by springs to the two particles that
        ###are before it and after it in the index list. This is always the case except
        ###for the last two and the first two particles.
        for i in range(self.floorSize, self.numParticles):
            #first particle in lattice no particles before it in the index
            if i == self.floorSize:
                ##adds springs between this particle and the next two particles
                ##other end particles are similar to this
                self.adjMatrix[i,[i+1,i+2]] = False
            #second particle in lattice
            elif i == self.floorSize + 1:
                self.adjMatrix[i,[i+1,i+2]] = False
                self.adjMatrix[i,i-1] = False
            #second to last particle in lattice
            elif i == self.numParticles-2:
                self.adjMatrix[i,i+1] = False
                self.adjMatrix[i,[i-1,i-2]] = False
            #last particle in lattice
            elif i == self.numParticles-1:
                self.adjMatrix[i,[i-1,i-2]] = False
            #middle particles in lattice
            else:
                self.adjMatrix[i,[i+1,i+2]] = False
                self.adjMatrix[i,[i-1,i-2]] = False
        return self.adjMatrix


    def summary(self):
        print self.xpos
        print self.ypos
        print self.zpos
        print self.xvel
        print self.yvel
        print self.zvel
        print self.xacl
        print self.yacl
        print self.zacl
        print self.massVector
        print self.numParticles
        print self.Lx
        print self.Ly
        print self.Lz

    def posSummary(self):
        print self.xpos
        print self.ypos
        print self.zpos
    ##reset container to initial conditions of all particles, useful for multiple
    ##simulations on one container
    def reset(self):
        self.xpos = self.xposinit
        self.ypos = self.yposinit
        self.zpos = self.zposinit
        self.xvel = self.xvelinit
        self.yvel = self.yvelinit
        self.zvel = self.zvelinit
        self.xacl = self.xaclinit
        self.yacl = self.yaclinit
        self.zacl = self.zaclinit
        self.dragForce = []
    
#######Initialize container from a file#########
def InitalizeContainer(self,filename):
    with open(filename, 'r') as f:
      reader = csv.reader(f,delimiter=" ")
      self.numParticles = int(reader.next()[0])
      xdim = float(reader.next()[0])
      ydim = float(reader.next()[0])
      zdim = float(reader.next()[0])
      container = Container(xdim,ydim,zdim)
      for row in reader:
          state = [float(i) for i in row]
          container.addParticle(state[0:-1],state[-1])
      return container
