# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 16:35:26 2013

@author: Nathan Sponberg
"""

from numpy import array, copy
from DistMatrix import *

class Integrate(object):
    def __init__(self):
        pass

class Verlet(Integrate):
    def __init__(self, dt):
        self.dt = dt
        pass

    def __call__(self, function, container):
        c = container
        f = function
        boundryCalc = DistMatrix()
        #boundryCalc.UpdateNeighbors(c, 2**(1/6))
        '''
        boundryCalc.UpdatePBoundries(c.xpos, c.Lx)
        boundryCalc.UpdatePBoundries(c.ypos, c.Ly)
        boundryCalc.UpdatePBoundries(c.zpos, c.Lz)
        '''
        c.particleFlux.append(0)
        c.xpos = c.xpos + c.xvel *self.dt + .5 * c.xacl * self.dt**2
        c.ypos = c.ypos + c.yvel *self.dt + .5 * c.yacl * self.dt**2
        c.zpos = c.zpos + c.zvel *self.dt + .5 * c.zacl * self.dt**2
        currentXacl = copy(c.xacl)
        currentYacl = copy(c.yacl)
        currentZacl = copy(c.zacl)
        nextXacl, nextYacl, nextZacl = f(c)
        c.xvel = c.xvel +.5*(nextXacl + currentXacl)*self.dt
        c.yvel = c.yvel +.5*(nextYacl + currentYacl)*self.dt
        c.zvel = c.zvel +.5*(nextZacl + currentZacl)*self.dt
        c.xacl = nextXacl
        c.yacl = nextYacl
        c.zacl = nextZacl
        c.integrationIteration += 1
        #########################################

#      c.x = c.x + c.vx * self.__dt + .5 * c.ax * self.__dt ** 2
