#!/usr/bin/env python

import AOPSO_singleRun as AOPSO
import pylab
from numpy import array, linspace, log10
import ipdb as pdb

def TestOverStep():
    nDim = 10
    numOfParticles = 20
    maxIteration = 1000
    minX = array([-6.0]*nDim)
    maxX = array([6.0]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    numOfTrial = 20
    intDim = 5
    for i in xrange(11):
        step = i*1
        gBest = array([0.0]*maxIteration)
        for i in xrange(numOfTrial):
            p1 = AOPSO.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, AOPSO.Rastrigrin,intDim,step)
            p1.run()
            gBest = gBest + p1.gBestArray[:maxIteration]
        gBest = gBest / numOfTrial
        pylab.plot(range(maxIteration), gBest,label='step='+str(step))
    pylab.title('Rastrigrin function, intDim=5 ($G_{best}$ over 20 trials)')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('Average gBest over '+str(numOfTrial)+' runs')
    pylab.grid(True)
    pylab.yscale('log')
    ylim = [-6, 1]
    ystep = 1.0
#    pylab.ylim(ylim[0], ylim[1])
#    yticks = linspace(ylim[0], ylim[1], int((ylim[1]-ylim[0])/ystep+1))
#    pylab.yticks(tuple(yticks), tuple(map(str,yticks)))
    pylab.legend(loc='upper right')
    pylab.show()

if __name__=='__main__':
    TestOverStep()
