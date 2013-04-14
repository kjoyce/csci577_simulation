# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 16:51:28 2013

@author: Nathan Sponberg
"""

from DistMatrix import *
import numpy as ny
from numpy import ma, triu, array, zeros, isnan
from pylab import norm, dot
from IPython.core.debugger import Tracer
debug_here = Tracer()

class Force(object):
    def __init__(self):
        pass

class LennardJonesForce(object):
    def __init__(self, sigma = 1., epsilon = 1.):
        self.sigma = sigma
        self.epsilon = epsilon
        pass

    def __call__(self, container):
        self.c = container
        distCalc = Neighbors()
        distx, disty, distz = distCalc.CalcDist(self.c.xpos,
                                                 self.c.ypos,
                                                 self.c.zpos)
        distx = distCalc.CheckDist(distx, self.c.Lx)
        disty = distCalc.CheckDist(disty, self.c.Ly)
        distz = distCalc.CheckDist(distz, self.c.Lz)
        #debug_here()
        distCalc.UpdateNeighbors(self.c)
        masses = self.c.massVector
        distr = ma.masked_array(sqrt(distx**2 + disty**2 + distz**2),
                                [distx**2 + disty**2 + distz**2 == 0])
        K1 = (ma.divide(self.sigma,distr).filled(0.))**12
        K2 = (ma.divide(self.sigma,distr).filled(0.))**6
        #debug_here()
        KE = sum(sum(triu(array(4.*self.epsilon*(K1-K2)))))
        print KE
        K3 = 2*K1 - K2
        magnitude = 24.*ma.divide(self.epsilon,distr)*K3
        xacl = ny.sum(array((magnitude * ma.divide(distx,distr)).filled(0.)), axis = 1)/masses
        yacl = ny.sum(array((magnitude * ma.divide(disty,distr)).filled(0.)), axis = 1)/masses
        zacl = ny.sum(array((magnitude * ma.divide(distz,distr)).filled(0.)), axis = 1)/masses
        return xacl, yacl, zacl
      
class LennardJonesNeighbors(object):
    def __init__(self, sigma = 1., epsilon = 1.):
        self.sigma = sigma
        self.epsilon = epsilon
        pass    
    
    def __call__(self,container):
        self.c = container
        distCalc = Neighbors()
        xacl = zeros(self.c.numParticles)
        yacl = zeros(self.c.numParticles)
        zacl = zeros(self.c.numParticles)
        masses = self.c.massVector
        print "tick"
        for particle in range(self.c.numParticles):
            neighbors = self.c.neighborList[particle]
            distx, disty, distz, distr = distCalc.NeighborDist(particle, 
                                                               neighbors, 
                                                               self.c)
            K1 = (ma.divide(self.sigma,distr).filled(0.))**12
            K2 = (ma.divide(self.sigma,distr).filled(0.))**6
            K3 = 2*K1 - K2
            magnitude = 24.*ma.divide(self.epsilon,distr)*K3
            #debug_here()
            xacl[particle] = sum(array((magnitude * ma.divide(distx,distr)).filled(0.)))/masses[particle]
            yacl[particle] = sum(array((magnitude * ma.divide(disty,distr)).filled(0.)))/masses[particle]
            zacl[particle] = sum(array((magnitude * ma.divide(distz,distr)).filled(0.)))/masses[particle]
            #debug_here()
        return xacl, yacl, zacl
        
class GranularForces(object):
    def __init__(self, sigma = 1., epsilon = 1., gamma = 30.):
        self.sigma = sigma
        self.epsilon = epsilon
        self.gamma = gamma
        pass    
    
    def __call__(self, container):
        self.c = container
        distCalc = Neighbors()
        #distCalc.UpdateNeighbors(self.c,2**(1/6))
        xacl = zeros(self.c.numParticles)
        yacl = zeros(self.c.numParticles)
        zacl = zeros(self.c.numParticles)
        masses = self.c.massVector
        #print self.c.neighborList[-1]
        for particle in range(self.c.numParticles):
            if self.c.ypos[particle] <= self.c.openingPosition and self.c.ypos[particle] >= (self.c.openingPosition - 2):
                self.c.particleFlux[self.c.integrationIteration] += 1
            neighbors = self.c.neighborList[particle]
            distx, disty, distz, distr, relVelx, relVely, relVelz = distCalc.NeighborDist(particle, 
                                                               neighbors, 
                                                               self.c)                
            K1 = (ma.divide(self.sigma,distr).filled(0.))**12
            K2 = (ma.divide(self.sigma,distr).filled(0.))**6
            K3 = 2*K1 - K2
            magnitude = 24*ma.divide(self.epsilon,distr)*K3
            #debug_here()
            #rUnitVector = distr/norm()
            dampingForcex = zeros(len(neighbors))
            dampingForcey = zeros(len(neighbors))
            dampingForcez = zeros(len(neighbors))            
            for i in range(len(neighbors)):
                
                    #print relVelx, relVely
                   #debug_here()
                displacement = array([distx[i],disty[i],distz[i]])
                unitVector = displacement/norm(displacement)
                nans = isnan(unitVector)
                unitVector[nans] = 0.
                dotProduct = dot(displacement, (relVelx[i],relVely[i],relVelz[i]))
                #debug_here()
                forceVector = -1*self.gamma*dotProduct*unitVector
                dampingForcex[i] = forceVector[0]
                dampingForcey[i] = forceVector[1]
                dampingForcez[i] = forceVector[2]
                #if len(self.c.neighborList[-1]) > 1 and particle == (self.c.numParticles - 1):
                    #print displacement
                """
                dampingForcex[i] = (self.gamma*dotProduct*unitVector * ma.divide(distx[i],distr[i])).filled(0.)
                dampingForcey[i] = (self.gamma*dotProduct*unitVector + ma.divide(disty[i],distr[i])).filled(0.)
                dampingForcez[i] = (self.gamma*dotProduct*unitVector + ma.divide(distz[i],distr[i])).filled(0.)
                """
            #debug_here()
            distFromWall1, xWall1, yWall1 = distCalc.lineDist(-1*sqrt(3),-1,20.,
                                                              self.c.xpos[particle], self.c.ypos[particle])
            distFromWall2, xWall2, yWall2 = distCalc.lineDist(sqrt(3),-1,-7.,
                                                              self.c.xpos[particle], self.c.ypos[particle])
            distFromWall3, xWall3, yWall3 = distCalc.lineDist(0,1,0,
                                                              self.c.xpos[particle], self.c.ypos[particle])
            wall1AclX = 0.
            wall1AclY = 0.
            wall2AclX = 0.
            wall2AclY = 0.
            wall3AclX = 0.
            wall3AclY = 0.
            
            if distFromWall1 <= 2**(1/6) and self.c.ypos[particle] >=10:
                K1 = (self.sigma/distFromWall1)**12
                K2 = (self.sigma/distFromWall1)**6
                K3 = 2*K1 - K2
                wallMag1 = 24*(self.epsilon/distFromWall1)*K3
                displacement1 = array([xWall1,yWall1])
                unitVector = displacement1/(displacement1**2)
                nans = isnan(unitVector)
                unitVector[nans] = 0.
                dotProduct = dot(displacement1, (self.c.xvel[particle],self.c.yvel[particle]))
                #debug_here()
                forceVector1 = -1*self.gamma*dotProduct*unitVector
                wall1AclX = forceVector1[0]
                wall1AclX += wallMag1*(xWall1/distFromWall1)
                wall1AclY = forceVector1[1]
                wall1AclY += wallMag1*(yWall1/distFromWall1)
            if distFromWall2 <= 2**(1/6) and self.c.ypos[particle] >=10:
                K1 = (self.sigma/distFromWall2)**12
                K2 = (self.sigma/distFromWall2)**6
                K3 = 2*K1 - K2
                wallMag2 = 24*(self.epsilon/distFromWall2)*K3
                displacement2 = array([xWall2,yWall2])
                unitVector = displacement2/(displacement2**2)
                nans = isnan(unitVector)
                unitVector[nans] = 0.
                dotProduct = dot(displacement2, (self.c.xvel[particle],self.c.yvel[particle]))
                #debug_here()
                forceVector2 = -1*self.gamma*dotProduct*unitVector
                wall2AclX = forceVector2[0]
                wall2AclX += wallMag2*(xWall2/distFromWall2)
                wall2AclY = forceVector2[1]
                wall2AclY += wallMag2*(yWall2/distFromWall2)
            if distFromWall3 <= 2**(1/6):
                K1 = (self.sigma/distFromWall3)**12
                K2 = (self.sigma/distFromWall3)**6
                K3 = 2*K1 - K2
                wallMag3 = 24*(self.epsilon/distFromWall3)*K3
                displacement3 = array([xWall3,yWall3])
                unitVector = displacement3/(displacement3**2)
                nans = isnan(unitVector)
                unitVector[nans] = 0.
                dotProduct = dot(displacement3, (self.c.xvel[particle],self.c.yvel[particle]))
                #debug_here()
                forceVector3 = -1*self.gamma*dotProduct*unitVector
                wall3AclX = forceVector3[0]
                wall3AclX += wallMag3*(xWall3/distFromWall3)
                wall3AclY = forceVector3[1]
                wall3AclY += wallMag3*(yWall3/distFromWall3)
            
            xacl[particle] = sum(array((magnitude * ma.divide(distx,distr)).filled(0.)))/masses[particle]
            xacl[particle] += sum(dampingForcex)/masses[particle]
            xacl[particle] += wall1AclX
            xacl[particle] += wall2AclX
            xacl[particle] += wall3AclX
            yacl[particle] = sum(array((magnitude * ma.divide(disty,distr)).filled(0.)))/masses[particle]
            yacl[particle] += sum(dampingForcey)/masses[particle]
            yacl[particle] += wall1AclY
            yacl[particle] += wall2AclY
            yacl[particle] += wall3AclY
            zacl[particle] = sum(array((magnitude * ma.divide(distz,distr)).filled(0.)))/masses[particle]
            zacl[particle] += sum(dampingForcex)/masses[particle]
            yacl[particle] += -2.
            #debug_here()
             ###set floor particles to stationary states
        #xacl[0:(self.c.floorSize+2*self.c.wallSize+2*self.c.slantSize)] = 0.
        #yacl[0:(self.c.floorSize+2*self.c.wallSize+2*self.c.slantSize)] = 0.
        
        zacl *= 0.
        print self.c.particleFlux[self.c.integrationIteration]
        return xacl, yacl, zacl
