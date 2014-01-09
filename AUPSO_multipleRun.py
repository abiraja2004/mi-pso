#!/usr/bin/env python

import AUPSO_singleRun as AUPSO
import pylab
from numpy import array, linspace, log10
import ipdb as pdb

def SphereTestOverIntDim():
    nDim = 10
    numOfParticles = 20
    maxIteration = 200
    minX = array([-100.0]*nDim)
    maxX = array([100.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    alpha = 0.5
    for intDim in xrange(0,11,2):
        gBest = array([0.0]*maxIteration)
        for i in xrange(numOfTrial):
            p1 = AUPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, AUPSO.Sphere,intDim,alpha)
            p1.run()
            gBest = gBest + p1.gBestArray[:maxIteration]
        gBest = gBest / numOfTrial
        pylab.plot(range(maxIteration), log10(gBest),label='intDim='+str(intDim))
    pylab.title('$G_{best}$ over 20 trials'+' alpha='+str(alpha))
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.grid(True)
    pylab.yscale('log')
    ylim = [-6, 1]
    ystep = 1.0
#    pylab.ylim(ylim[0], ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])/ystep+1))
#    pylab.yticks(tuple(yticks), tuple(map(str,yticks)))
    pylab.legend()
    pylab.show()
   
def SphereTestOverAlpha():
    nDim = 10
    numOfParticles = 20
    maxIteration = 200
    minX = array([-100.0]*nDim)
    maxX = array([100.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    intDim = 5
    alpha = 1.0
    while alpha>1e-3:
        gBest = array([0.0]*maxIteration)
        for i in xrange(numOfTrial):
            p1 = AUPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, AUPSO.Sphere,intDim,alpha)
            p1.run()
            gBest = gBest + p1.gBestArray[:maxIteration]
        gBest = gBest / numOfTrial
        pylab.plot(range(maxIteration), log10(gBest),label='alpha='+str(alpha))
        print 'alpha = ', alpha
        alpha -= 0.2
    print 'now drawing'
    pylab.title('$G_{best}$ over 20 trials'+' nDim='+str(nDim))
    pylab.xlabel('The $N^{th}$ Iteratioin')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.grid(True)
    pylab.yscale('log')
    ylim = [-6, 1]
    ystep = 1.0
#    pylab.ylim(ylim[0], ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])/ystep+1))
#    pylab.yticks(tuple(yticks), tuple(map(str,yticks)))
    pylab.legend()
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
    pylab.legend()
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
    SphereTestOverAlpha()
