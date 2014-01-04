#!/usr/bin/env python

import RPSO_singleRun as RPSO
import pylab
from numpy import array, linspace
import ipdb as pdb

def SphereTest():
    nDim = 10
    numOfParticles = 5
    maxIteration = 500
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    gBest = array([0.0]*maxIteration)
    for i in xrange(numOfTrial):
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Sphere)
        p1.run()
        gBest = gBest + p1.gBestArray[:maxIteration]
    gBest = gBest / numOfTrial
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.grid(True)
    pylab.yscale('log')
#    ylim = [-0.5, 4.5]
#    ystep = 0.5
#    pylab.ylim(10**ylim[0], 10**ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int(ylim[1]-ylim[0]/ystep+1))
#    pylab.yticks(tuple(10**yticks), tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), gBest,'-')
    pylab.show()
   
        
def RosenbrockTest():
    nDim = 10
    numOfParticles = 5
    maxIteration = 500
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    gBest = array([0.0]*maxIteration)
    for i in xrange(numOfTrial):
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Rosenbrock)
        p1.run()
        gBest = gBest + p1.gBestArray[:maxIteration]
    gBest = gBest / numOfTrial
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
    pylab.yscale('log')
    pylab.ylim(10**-0.5, 10**4.5)
    yticks = linspace(-0.5,4.5, 5/0.5+1)
    pylab.yticks(tuple(10**yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), gBest,'-')
    pylab.show()
   
def RastrigrinTest():
    nDim = 3
    numOfParticles = 20
    maxIteration = 200
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    gBest = array([0.0]*maxIteration)
    for i in xrange(numOfTrial):
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Rastringrin)
        p1.run()
        gBest = gBest + p1.gBestArray[:maxIteration]
    gBest = gBest / numOfTrial
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
    pylab.yscale('log')
    pylab.ylim(10**-1.5, 10**1.5)
    yticks = linspace(-1.5,1.5, 3/0.5+1)
    pylab.yticks(tuple(10**yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), gBest,'-')
    pylab.show()

def GriewankTest():
    nDim = 3
    numOfParticles = 20
    maxIteration = 200
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    gBest = array([0.0]*maxIteration)
    for i in xrange(numOfTrial):
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Griewank)
        p1.run()
        gBest = gBest + p1.gBestArray[:maxIteration]
    gBest = gBest / numOfTrial
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
    pylab.yscale('log')
    pylab.ylim(10**-2.2, 10**-0.6)
    yticks = linspace(-2.2,-0.6, (2.2-0.6)/0.2+1)
    pylab.yticks(tuple(10**yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), gBest,'-')
    pylab.show()




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
    RosenbrockTest()
