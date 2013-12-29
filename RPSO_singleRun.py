#!/usr/bin/env python

from math import *
from random import random
import ipdb as pdb
from numpy import array
import pylab

class Particle:
    def __init__(self, nDim, minX=None, maxX=None, minV=None, maxV=None):
        self.nDim = nDim
        if minX== None:
            self.minX = [0]*nDim
        else:
            self.minX = minX
        if maxX==None:
            self.maxX = [0]*nDim
        else:
            self.maxX = maxX
        if minV==None:
            self.minV = [0]*nDim
        else:
            self.minV = minV
        if maxV==None:
            self.maxV = [0]*nDim
        else:
            self.maxV = maxV
        # v: velocity
        self.v = [0.0] * nDim
        # x: position value
        self.x = [0.0] * nDim
        for i in xrange(self.nDim):
            self.v[i] = random()*(self.maxV[i]-self.minV[i]) + self.minV[i]
            self.x[i] = random()*(self.maxX[i]-self.minX[i]) + self.minX[i]
        self.fit = self.fitness()
    # setPositiion not used in BPSO
    def setVelocity(self, v):
        self.v = v[:]
        for i in xrange(self.nDim):
            if self.v[i] < self.minV[i] :
                self.v[i] = self.minV[i]
            elif self.v[i] > self.maxV[i]:
                self.v[i] = self.maxV[i]
    def updatePosition(self, boundaryType='Invisible'):
        for i in xrange(self.nDim):
            self.x[i] += self.v[i]
            if (self.x[i] < self.minX[i]): #exceed up bound
                # three kinds of boundary handling
                if boundaryType=='Absorbing':
                    self.x[i] = self.minX[i]
                    self.fit = self.fitness()
                elif boundaryType=='Reflecting':
                    self.x[i] = self.minX[i] + self.minX[i]-self.x[i] 
                    self.fit = self.fitness()
                elif boundaryType=='Invisible':
                    self.fit = 1e50
            elif (self.x[i] > self.maxX[i]): #exceed down bound
                if boundaryType=='Absorbing':
                    self.x[i] = self.maxX[i]
                    self.fit = self.fitness()
                elif boundaryType=='Reflecting':
                    self.x[i] = self.maxX[i] + self.maxX[i]-self.x[i] 
                    self.fit = self.fitness()
                elif boundaryType=='Invisible':
                    self.fit = 1e50
            # position in boundary
            else:
                self.fit = self.fitness()
    def fitness(self):
        x = self.x[:]
        return Rastringrin(x)

class PSOProblem:
    # minX, maxX, minV, maxV should be list
    def __init__(self, nDim, numOfParticles=None, maxIteration=None, minX=None, maxX=None, minV = None, maxV = None):
        self.nDim = nDim
        if numOfParticles == None:
            self.numOfParticles = 10 
        else:
            self.numOfParticles = numOfParticles
        if maxIteration == None:
            self.maxIteration = 200
        else:
            self.maxIteration = maxIteration
        self.minX = minX
        self.maxX = maxX
        self.minV = minV
        self.maxV = maxV
        self.p = [None]*(self.numOfParticles)
        self.pBest = [None]*(self.numOfParticles)
        self.gBest = None
        self.initParticles()
    def initParticles(self):
        gBest = 1e50
        bestK = 0
        for k in xrange(self.numOfParticles):
            self.p[k] = Particle(self.nDim, self.minX,self.maxX,self.minV,self.maxV)
            self.pBest[k] = Particle(self.nDim, self.minX,self.maxX,self.minV,self.maxV)
            self.pBest[k].x = self.p[k].x[:]
            self.pBest[k].fit = self.p[k].fit
            if self.pBest[k].fit < gBest:
                gBest = self.pBest[k].fit
                bestK = k
        self.gBest = Particle(self.nDim)
        self.gBest.x = self.p[bestK].x[:]
        self.gBest.fit = gBest
#        print 'initial gBest fitness is %', self.gBest.fit
            
    def run(self):
        # w decreace from 0.9 to 0.4 over iterations
        w = 0.9
        c1 = 2.0
        c2 = 2.0
        # gBestArray: save gBest in each iteration
        # gBestArray[0]: initial gBest
        self.gBestArray = [None]*(self.maxIteration+1)
        self.gBestArray[0] = self.gBest.fit
        for i in xrange(self.maxIteration):
            for j in xrange(self.numOfParticles):
                w = 0.9-(0.5)/self.maxIteration*i
                v = [0.0]*(self.nDim)
                for dim in range(self.nDim):
                    v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                            + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                self.p[j].setVelocity(v)
                self.p[j].updatePosition()
                if isBetterThan(self.p[j], self.pBest[j]):
                    self.pBest[j].x = self.p[j].x[:]
                    self.pBest[j].fit = self.p[j].fit
                if isBetterThan(self.p[j], self.gBest):
                    self.gBest.x = self.p[j].x[:]
                    self.gBest.fit = self.p[j].fit
            self.gBestArray[i+1] = self.gBest.fit
    
    def drawResult(self):
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

# Test functions
def Sphere(x):
    result = 0.0
    for i in xrange(len(x)):
        result += x[i]**2 
    return result

def Rosenbrock(x):
    result = 0.0
    for i in xrange(len(x)-1):
        result += (1.0-x[i])**2 + 100.0*((x[i+1]-x[i]**2)**2)
    return result

def Rastringrin(x):
    result = 0.0
    for i in xrange(len(x)):
        result += x[i]**2 - 10.0*cos(2*pi*x[i])+10
    return result

def RPSOTest():
    nDim = 3
    numOfParticles = 10
    maxIteration = 200
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    p1 = PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV)
    p1.run()
    p1.drawResult()

if __name__=='__main__':
    RPSOTest()
