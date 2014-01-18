#!/usr/bin/env python

import RPSO_singleRun as RPSO
import BPSO_singleRun as BPSO
import pylab
from numpy import array, linspace, log10
import ipdb as pdb

def BPSOvsRPSOTest():
    nDim = 3
    numOfParticles = 10
    maxIteration = 200
    minX = array([-5.0]*nDim)
    maxX = array([5.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    lowerBound = -5.0
    upperBound =  5.0
    nSeg = 3
    BPSODim = 24
    gBest1 = array([0.0]*maxIteration)
    gBest2 = array([0.0]*maxIteration)
    numOfTrial = 200
    for i in xrange(numOfTrial):
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Rastrigrin)
        p1.run()
        gBest1 = gBest1 + p1.gBestArray[:maxIteration]
        p2 = BPSO.PSOProblem(BPSODim, numOfParticles, maxIteration,nSeg, lowerBound, upperBound, RPSO.Rastrigrin)
        p2.run()
        gBest2 = gBest2 + p2.gBestArray[:maxIteration]
    gBest1 = gBest1 / numOfTrial
    gBest2 = gBest2 / numOfTrial
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.title('$G_{best}$ over 20 trials')
    pylab.grid(True)
#    pylab.yscale('log')
#    ylim = [-6, 1]
    ystep = 1.0
#    pylab.ylim(ylim[0], ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])/ystep+1))
#    pylab.yticks(tuple(yticks), tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), gBest1, label='RPSO')
    pylab.plot(range(maxIteration), gBest2, label='BPSO')
    pylab.legend()
    pylab.show()
   
def TestOverMaxV():
    nDim = 10
    numOfParticles = 20
    maxIteration = 200
    minX = array([-100.0]*nDim)
    maxX = array([100.0]*nDim)
    numOfTrial = 20
    coef = 0.1
    while coef<1.0:
        maxV = coef*(maxX - minX)
        minV = -1.0*maxV
        gBest = array([0.0]*maxIteration)
        for i in xrange(numOfTrial):
            p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Sphere)
            p1.run()
            gBest = gBest + p1.gBestArray[:maxIteration]
        gBest = gBest / numOfTrial
        pylab.plot(range(maxIteration), log10(gBest),label='coef='+str(coef))
        coef = coef+0.2
    pylab.title('$G_{best}$ over 20 trials')
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.grid(True)
#    pylab.yscale('log')
    ylim = [-6, 1]
    ystep = 1.0
#    pylab.ylim(ylim[0], ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])/ystep+1))
#    pylab.yticks(tuple(yticks), tuple(map(str,yticks)))
    pylab.legend(loc='lower left')
    pylab.show()
        

def RosenbrockTest():
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
        p1 = RPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, RPSO.Rosenbrock)
        p1.run()
        gBest = gBest + p1.gBestArray[:maxIteration]
    gBest = gBest / numOfTrial
    pylab.title('$G_{best}$ over 20 trials')
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
#    pylab.yscale('log')
    ymin, ymax = -1.5, 2.5
    ystep = 0.5
    pylab.ylim(ymin, ymax)
    yticks = linspace(ymin, ymax, (ymax-ymin)/ystep+1)
    pylab.yticks(tuple(yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), log10(gBest),'-', label='Global best')
    pylab.legend(loc='upper right')
    pylab.show()
   
def RastringrinTest():
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
    pylab.title('$G_{best}$ over 20 trials')
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
#    pylab.yscale('log')
    ymin, ymax = -1.5, 1.5
    ystep = 0.5
    pylab.ylim(ymin, ymax)
    yticks = linspace(ymin, ymax, (ymax-ymin)/ystep+1)
    pylab.yticks(tuple(yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), log10(gBest),'-', label='Global best')
    pylab.legend()
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
    pylab.title('$G_{best}$ over 20 trials')
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs (logscale)')
    pylab.grid(True)
#    pylab.yscale('log')
    ymin, ymax = -2.2, -0.6
    ystep = 0.2
    pylab.ylim(ymin, ymax)
    yticks = linspace(ymin, ymax, (ymax-ymin)/ystep+1)
    pylab.yticks(tuple(yticks),tuple(map(str,yticks)))
    pylab.plot(range(maxIteration), log10(gBest),'-', label='Global best')
    pylab.legend()
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
    BPSOvsRPSOTest()
