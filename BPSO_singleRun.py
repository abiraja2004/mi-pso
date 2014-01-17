#!/usr/bin/env python

from testfunc import *
from math import *
from random import random
from numpy import array
from copy import deepcopy
from testfunc import *
import ipdb as pdb
import pylab

class Particle:
    """
    nDim: length of binary string
    nSeg: binary string are seperated into nSeg segment
    lowerBound and upperBound: each segment are transformed to [-lowerBound, upperBound]
    x: binary string
    """
    def __init__(self, nDim, nSeg,lowerBound, upperBound, fitfunc=None):
        self.nDim = nDim
        self.nSeg = nSeg
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.fitfunc = fitfunc
        # v: velocity
        self.v = [0.0] * nDim
        # x: binary value
        self.x = [0] * nDim
        for i in range(self.nDim):
            #According to Nanbo Jin's paper, v_max = 6
            self.v[i] = 0.0
        self.updateX()
        self.fit = self.fitness()
    def updateX(self):
        for i in range(self.nDim):
            s = 1.0 / (1.0 + exp(-self.v[i]))
            prob = random()
            if (prob < s) :
                self.x[i] = 1.0
            else:
                self.x[i] = 0.0
    # setPositiion not used in BPSO
    def setVelocity(self, v):
        self.v = v[:]
        for i in range(self.nDim):
            if self.v[i] > 6 :
                self.v[i] = 6.0
            elif self.v[i] < -6:
                self.v[i] = -6.0
    def updatePosition(self):
        self.updateX()
        self.fit = self.fitness()
    def fitness(self):
        realX = [0.0] * self.nSeg
        segLen = len(self.x)/self.nSeg 
        coef = 1.0/(pow(2,segLen)-1)*(self.upperBound - self.lowerBound)
        for i in xrange(self.nSeg):
            bitStr = ''.join(map(str,map(int,(self.x[i*segLen:(i+1)*segLen]))))
            realX[i] = coef*int(bitStr,2)+self.lowerBound
        return self.fitfunc(realX)

class PSOProblem:
    def __init__(self, nDim, numOfParticles, maxIteration, nSeg, lowerBound, upperBound, fitfunc):
        self.nDim = nDim
        if numOfParticles == None:
            self.numOfParticles = 10 
        else:
            self.numOfParticles = numOfParticles
        if maxIteration == None:
            self.maxIteration = 200
        else:
            self.maxIteration = maxIteration
        self.nSeg = nSeg
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.fitfunc = fitfunc
        self.p = [0]*(self.numOfParticles)
        self.pBest = [0]*(self.numOfParticles)
        self.gBest = None
        self.initParticles()

    def initParticles(self):
        gBest = 1e50
        bestK = 0
        for k in range(self.numOfParticles):
            self.p[k] = Particle(self.nDim, self.nSeg, self.lowerBound, self.upperBound, self.fitfunc)
            self.pBest[k] = deepcopy(self.p[k])
            self.pBest[k].x = self.p[k].x[:] #This is not necessary
            self.pBest[k].fit = self.p[k].fit
            if self.pBest[k].fit < gBest:
                gBest = self.pBest[k].fit
                bestK = k
        self.gBest = deepcopy(self.p[bestK])
        self.gBest.x = self.p[bestK].x[:]
        self.gBest.fit = gBest
#        print 'initial gBest fitness is %', self.gBest.fit
            
    def run(self):
        w = 1.0
        c1 = 2.0
        c2 = 2.0
        self.gBestArray = array([0.0]*(self.maxIteration+1))
        self.gBestArray[0] = self.gBest.fit
        for i in range(self.maxIteration):
            for j in range(self.numOfParticles):
                #pdb.set_trace()
                v = [0]*(self.nDim)
                for dim in range(self.nDim):
                    v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                            + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                self.p[j].setVelocity(v)
#                print 'in run, v', self.p[j].v
                self.p[j].updatePosition()
#                print 'in run', self.p[j].fit
                if isBetterThan(self.p[j], self.pBest[j]):
                    self.pBest[j].x = self.p[j].x[:]
                    self.pBest[j].fit = self.p[j].fit
                if isBetterThan(self.p[j], self.gBest):
                    self.gBest.x = self.p[j].x[:]
                    self.gBest.fit = self.p[j].fit
            self.gBestArray[i+1] = self.gBest.fit

    def drawResult(self):
        print 'Init gBest:', self.gBestArray[0]
        print 'Final gBest:', self.gBestArray[self.maxIteration]
        pylab.plot(range(1+self.maxIteration), self.gBestArray,'-')
        pylab.xlabel('The $N^{th}$ Iteratioin')
        pylab.ylabel('Global Best')
        pylab.grid(True)
        pylab.show()


# smaller is better
def isBetterThan(particleA, particleB):
    result = None
    if (particleA.fit < particleB.fit) :
        result = True
    else:
        result = False
    return result

def RastrigrinTest():
    nDim = 48
    nSeg = 3
    lowerBound = -5.0
    upperBound = 5.0
    numOfParticles = 10
    maxIteration = 200
    p1 = PSOProblem(nDim, numOfParticles, maxIteration, nSeg, lowerBound, upperBound, Rastrigrin)
    p1.run()
    p1.drawResult()

if __name__=='__main__':
    RastrigrinTest()
