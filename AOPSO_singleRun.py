#!/usr/bin/env python

from testfunc import *
from math import *
from random import random
import ipdb as pdb
from numpy import array
import pylab

# x[0:-intDim]: real parameter
# x[-intDim:]: integer parameter
class Particle:
    def __init__(self, nDim, minX=None, maxX=None, minV=None, maxV=None, fitfunc=None, penaltyfunc=None, intDim = 0):
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
        self.fitfunc = fitfunc
        self.penaltyfunc = penaltyfunc
        self.intDim = intDim
        # Store current iteration number in optimization
        self.iternum= 0
        # v: velocity
        self.v = [0.0] * nDim
        # x: position value
        self.x = [0.0] * nDim
        for i in xrange(self.nDim):
            self.v[i] = random()*(self.maxV[i]-self.minV[i]) + self.minV[i]
            self.x[i] = random()*(self.maxX[i]-self.minX[i]) + self.minX[i]
        for i in xrange(self.nDim-self.intDim, self.nDim):
            self.x[i] = round(self.x[i])
        self.fit = self.fitness()
    def setIternum(self, k):
        self.iternum = k
    def getIternum(self):
        return self.iternum
    # setPositiion not used in BPSO
    def setVelocity(self, v):
        self.v = v[:]
        for i in xrange(self.nDim):
            if (i>=self.nDim-self.intDim): # discrete variable's velocity is also integer
                self.v[i] = round(self.v[i])
            if self.v[i] < self.minV[i] :
                self.v[i] = self.minV[i]
            elif self.v[i] > self.maxV[i]:
                self.v[i] = self.maxV[i]

    def updatePosition(self, boundaryType='Invisible'):
        worstFitness = False
        for i in xrange(self.nDim):
            if (i < (self.nDim-self.intDim)):
                self.x[i] += self.v[i]
            else:
                self.x[i] += self.v[i]
                self.x[i] = round(self.x[i]) # Enforce round on integer parameter
            if (self.x[i] < self.minX[i]): #exceed up bound
                # three kinds of boundary handling
                if boundaryType=='Absorbing':
                    self.x[i] = self.minX[i]
                elif boundaryType=='Reflecting':
                    self.x[i] = self.minX[i] + self.minX[i]-self.x[i] 
                elif boundaryType=='Invisible':
                    worstFitness = True
            elif (self.x[i] > self.maxX[i]): #exceed down bound
                if boundaryType=='Absorbing':
                    self.x[i] = self.maxX[i]
                elif boundaryType=='Reflecting':
                    self.x[i] = self.maxX[i] + self.maxX[i]-self.x[i] 
                elif boundaryType=='Invisible':
                    worstFitness = True
        if (worstFitness):
            self.fit = 1e50
        else:
            self.fit = self.fitness()
    def fitness(self):
        x = self.x[:]
        fit = self.fitfunc(x)
        if (self.penaltyfunc==None):
            penalty = 0.0
        else:
            penalty = self.penaltyfunc(self.iternum,x)
        return fit + penalty

class PSOProblem:
    """
    Alternating Optimized PSO
    continuous and discrete variables are updated alternatively
    step = 0: Just the same as ordinary PSO
    step = 1: continuous, discrete, continuous, discrete,...
    stopstep: after stopstep, stop alternating
    """
    # minX, maxX, minV, maxV should be list/array
    def __init__(self, nDim, numOfParticles=None, maxIteration=None, minX=None, maxX=None, minV = None, maxV = None, fitfunc=None, penaltyfunc=None, intDim=0, step=0, stopstep=1e50 ):
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
        self.fitfunc = fitfunc
        self.penaltyfunc = penaltyfunc
        self.intDim = intDim
        self.step = step
        self.stopstep = stopstep
        self.p = [None]*(self.numOfParticles)
        self.pBest = [None]*(self.numOfParticles)
        self.gBest = None
        self.initParticles()
    def initParticles(self):
        gBest = 1e50
        bestK = 0
        for k in xrange(self.numOfParticles):
            self.p[k] = Particle(self.nDim, self.minX,self.maxX,self.minV,self.maxV,self.fitfunc,self.penaltyfunc, self.intDim)
            self.pBest[k] = Particle(self.nDim, self.minX,self.maxX,self.minV,self.maxV,self.fitfunc, self.penaltyfunc, self.intDim)
            self.pBest[k].x = self.p[k].x[:]
            self.pBest[k].fit = self.p[k].fit
            if self.pBest[k].fit < gBest:
                gBest = self.pBest[k].fit
                bestK = k
        self.gBest = Particle(self.nDim, fitfunc=self.fitfunc)
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
        self.gBestArray = array([0.0]*(self.maxIteration+1))
        self.gBestArray[0] = self.gBest.fit
        updateIntDim = True
        for i in xrange(self.maxIteration):
            for j in xrange(self.numOfParticles):
                w = 0.9-(0.5)/self.maxIteration*i
                v = [0.0]*(self.nDim)
                if (self.step==0 or i>self.stopstep):
                    for dim in range(self.nDim):
                        v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                                + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                else:
                    if (updateIntDim):
                        for dim in range(self.nDim-self.intDim, self.nDim):
                            v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                                    + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                    else:
                        for dim in range(self.nDim-self.intDim):
                            v[dim] = ( w * self.p[j].v[dim] + c1*random()*(self.pBest[j].x[dim] - self.p[j].x[dim])
                                    + c2*random()*(self.gBest.x[dim] - self.p[j].x[dim]))
                    if (i%self.step==0):
                        updateIntDim = not updateIntDim
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
        pylab.title('PSO')
        pylab.xlabel('The $N^{th}$ Iteratioin')
        pylab.ylabel('Global Best')
        pylab.yscale('log')
        pylab.grid(True)
        pylab.plot(range(1+self.maxIteration), self.gBestArray,'-', label=self.fitfunc.__name__)
        pylab.legend()
        pylab.show()


# smaller is better
def isBetterThan(particleA, particleB):
    result = None
    if (particleA.fit < particleB.fit) :
        result = True
    else:
        result = False
    return result


#nD sinc 
def Sinc(x):
    result = 1.0
    myerr = 1e-8
    for i in xrange(len(x)):
        if (x[i] > myerr):
            result *= sin(pi*x[i]) / (pi*x[i])
        else:
            result *= sin(pi*myerr)/ (pi*myerr)
    result = 1.0-result
    return result

def SphereTest():
    nDim = 10
    numOfParticles = 20
    maxIteration = 2000
    minX = array([-100.0]*nDim)
    maxX = array([100.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    step = 0
    p = [None]*5
    for i in xrange(5):
        step = i*5
        p[i] = PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, Sphere, 5, step)
        p[i].run()
    pylab.title('AO-PSO')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('Global Best')
    pylab.grid(True)
    pylab.yscale('log')
    for i in xrange(5):
        pylab.plot(range(1+maxIteration), p[i].gBestArray,'-', label='step='+str(i*5))
    pylab.legend(loc='lower left')
    pylab.show()

def RosenbrockTest():
    nDim = 10
    numOfParticles = 20
    maxIteration = 2000
    minX = array([-10.0]*nDim)
    maxX = array([10.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    step = 0
    p = [None]*5
    for i in xrange(5):
        step = i*5
        p[i] = PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, Rosenbrock, intDim=5, step=step)
        p[i].run()
    pylab.title('AO-PSO Rosenbrock')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('Global Best')
    pylab.grid(True)
    pylab.yscale('log')
    for i in xrange(5):
        pylab.plot(range(1+maxIteration), p[i].gBestArray,'-', label='step='+str(i*5))
    pylab.legend(loc='lower left')
    pylab.show()





def GriewankTest():
    nDim = 10
    intDim = 3
    numOfParticles = 10
    maxIteration = 200
    minX = array([-600.0]*nDim)
    maxX = array([600.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    alpha = 0.5
    p1 = PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV,Griewank, intDim, alpha)
    p1.run()
    alpha = 0.0
    p2 = PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV,Griewank, intDim, alpha)
    p2.run()
    pylab.title('AU-PSO')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('Global Best')
    pylab.grid(True)
    pylab.plot(range(1+maxIteration), p1.gBestArray,'-', label='alpha=0.5')
    pylab.plot(range(1+maxIteration), p2.gBestArray,'-', label='alpha=0.0')
    pylab.legend(loc='upper right')
    pylab.show()




if __name__=='__main__':
    RosenbrockTest()
