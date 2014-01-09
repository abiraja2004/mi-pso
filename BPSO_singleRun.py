#!/usr/bin/env python

from testfunc import *
from math import *
from random import random
import ipdb as pdb

class Particle:
    def __init__(self, nDim):
        self.nDim = nDim
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
        x = self.x[:]
        nDim = self.nDim
        a = [0.0] * 3
        fitness = 0.0
        for i in range(16):
            a[0] = a[0] + 1.0*x[i] * pow(2,i)
        for i in range(16,32):
            a[1] = a[1] + 1.0*x[i] * pow(2,i-16)
        for i in range(32, 48):
            a[2] = a[2] + 1.0*x[i] * pow(2,i-32)
        for i in range(3):
            a[i] = a[i] / (pow(2,16)-1.0) * 10.0 - 5.0
            fitness = fitness + a[i]*a[i] - 10*cos(2*pi*a[i]) + 10.0
        return fitness

class PSOProblem:
    # integerDim: store dims in which variables should be integers.
    def __init__(self, nDim, numOfParticles=None, maxIteration=None):
        self.nDim = nDim
        if numOfParticles == None:
            self.numOfParticles = 10 
        else:
            self.numOfParticles = numOfParticles
        if maxIteration == None:
            self.maxIteration = 200
        else:
            self.maxIteration = maxIteration
        self.p = [0]*(self.numOfParticles)
        self.pBest = [0]*(self.numOfParticles)
        self.gBest = None
        self.initParticles()

    def initParticles(self):
        gBest = 1e50
        bestK = 0
        for k in range(self.numOfParticles):
            self.p[k] = Particle(self.nDim)
            self.pBest[k] = Particle(self.nDim)
            self.pBest[k].x = self.p[k].x[:]
        #    self.pBest[k].v = self.p[k].v[:]
            self.pBest[k].fit = self.p[k].fit
            if self.pBest[k].fit < gBest:
                gBest = self.pBest[k].fit
                bestK = k
        self.gBest = Particle(self.nDim)
        self.gBest.x = self.p[bestK].x[:]
        self.gBest.fit = gBest
#        print 'initial gBest fitness is %', self.gBest.fit
            
    def run(self):
        self.gStat = [0] * 200
        w = 1.0
        c1 = 2.0
        c2 = 2.0
        self.initParticles()
        print 'Iteration \tgBest \tR1 \tR2 \tR3'
        for i in range(self.maxIteration):
            for j in range(self.numOfParticles):
                # update velocity
                #pdb.set_trace()
                v = [0]*(self.nDim)
                for dim in range(self.nDim):
                    v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                            + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                self.p[j].setVelocity(v)
#                print 'in run, v', self.p[j].v
                self.p[j].updatePosition()
#                print 'in run, x', self.p[j].x
#                print 'in run', self.p[j].fit
                if isBetterThan(self.p[j], self.pBest[j]):
                    self.pBest[j].x = self.p[j].x[:]
                    self.pBest[j].fit = self.p[j].fit
                if isBetterThan(self.p[j], self.gBest):
                    self.gBest.x = self.p[j].x[:]
                    self.gBest.fit = self.p[j].fit
            #pdb.set_trace()
            print i, '\t', self.gBest.fit, '\t', getR1(self.gBest.x), '\t', getR2(self.gBest.x), '\t', getR3(self.gBest.x)
#            print 'Round', i, 'gbest', self.gBest.fit
#                print i, self.gBest.fit

# smaller is better
def isBetterThan(particleA, particleB):
    result = None
    if (particleA.fit < particleB.fit) :
        result = True
    else:
        result = False
    return result

# convert 16 bit binary number to real in [-5, 5]
def binaryToReal(x):
    result = 0
    for i in range(16):
        result += pow(2,i)*x[i]
    result = result / (pow(2,16)-1)* 10 - 5.0
    return result

# R1: x[47]..x[32]
# R2: x[31]..x[16]
# R3: x[15]..x[0]
def getR1(x):
    return ''.join(map(str,map(int,x[32:48][::-1])))+'(%.3g)'%(binaryToReal(x[32:48]))

def getR2(x):
    return ''.join(map(str,map(int,x[16:32][::-1])))+'(%.3g)'%(binaryToReal(x[16:32]))

def getR3(x):
    return ''.join(map(str,map(int,x[:16][::-1])))+'(%.3g)'%(binaryToReal(x[:16]))

def BPSOTest():
    nDim = 48
    numOfParticles = 80
    maxIteration = 200
    p1 = PSOProblem(nDim, numOfParticles, maxIteration)
    p1.run()

if __name__=='__main__':
    BPSOTest()
