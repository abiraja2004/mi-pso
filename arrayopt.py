#!/usr/bin/env python

import AOPSO_singleRun as aopso
import BPSO_singleRun as bpso
import RPSO_singleRun as rpso
from arrayfactor import OneDimAF, OneDimSLL
from numpy import linspace, abs, sum, pi, sin, exp, log10, max, zeros, array,concatenate
import pylab
import ipdb as pdb

def use_aopso():
    nDim = 10
    numOfParticles = 20
    maxIteration = 1000
    minX = array([0.125]*nDim)
    maxX = array([4.625]*nDim)
    minX[9] = 5
    maxX[9] = 511
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    maxV[9] = int(maxV[9])
    minV[9] = int(minV[9])
    step=0
    p1 = aopso.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, arrayfit, arraypenalty, intDim=1, step=step)
    p1.run()
    pylab.title('array optimization')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('SLL')
    pylab.grid(True)
    pylab.ylim(-25,0)
    pylab.plot(range(1+maxIteration), p1.gBestArray, label='SLL')
    pylab.legend()
    pylab.show()
    #output array config
    print getArrayData(p1.gBest.x)
    
def use_rpso():
    nDim = 8
    numOfParticles = 20
    maxIteration = 1000
    minX = array([0.125]*nDim)
    maxX = array([4.625]*nDim)
    maxV = 0.2*(maxX - minX)
    minV = -1.0*maxV
    p1 = rpso.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, arrayfit2, arraypenalty2)
    p1.run()
    pylab.title('array optimization')
    pylab.xlabel('The $N^{th}$ Iteration')
    pylab.ylabel('SLL')
    pylab.grid(True)
    pylab.ylim(-25,0)
    pylab.plot(range(1+maxIteration), p1.gBestArray, label='SLL')
    pylab.legend()
    pylab.show()
    #output array config
    arrayData = p1.gBest.x[:]
    arrayData.sort()
    print arrayData




def getArrayData(x):
    elemNum = 9
    state = getBin(int(x[elemNum]),elemNum)
    arrayData = []
    for i in xrange(elemNum):
        if state[i]=='1':
            arrayData.append(x[i])
    arrayData = array(arrayData)
    arrayData.sort()
    return arrayData


def arraypenalty(k, x):
    """Return penalty value for x in kth iteration
    Penalty = 1.0* numOfCloseInterval
    """
    elemNum = 9
    state = getBin(int(x[elemNum]),elemNum)
    arrayData = []
    for i in xrange(elemNum):
        if state[i]=='1':
            arrayData.append(x[i])
    arrayData = array(arrayData)
    arrayData.sort()
    numOfInvalidInterval = 0
    if (len(arrayData)>1):
        for i in xrange(len(arrayData)-1):
            if (arrayData[i+1] - arrayData[i])<0.25:
                numOfInvalidInterval+=1
    return 10.0*numOfInvalidInterval

def arrayfit(x):
    """
    x[0:9]: position of array elements
    x[9]:  on/off state of 9 array elements (0-511)
    """
    elemNum = 9
    state = getBin(int(x[elemNum]), elemNum)
    arrayData = []
    for i in xrange(elemNum):
        if state[i]=='1':
            arrayData.append(x[i])
    arrayData.append(4.75)
    arrayData = array(arrayData)
    pos = concatenate((-1*(arrayData[::-1]),arrayData))
    num = len(pos)
    amp = 1.0
    phase = 0.0
    AF = zeros(360)
    for i in xrange(360):
        AF[i] = OneDimAF(num, pos, amp, phase,i/180.*pi) 
    AF = AF - max(AF)
    return OneDimSLL(AF)

def arrayfit2(x):
    # optimize 8 elem
    arrayData = x[:]
    arrayData.append(4.75)
    arrayData = array(arrayData)
    pos = concatenate((-1*(arrayData[::-1]),arrayData))
    num = len(pos)
    amp = 1.0
    phase = 0.0
    AF = zeros(360)
    for i in xrange(360):
        AF[i] = OneDimAF(num, pos, amp, phase,i/180.*pi) 
    AF = AF - max(AF)
    return OneDimSLL(AF)

def arraypenalty2(k, x):
    """Return penalty value for 8 elem in kth iteration
    Penalty = 1.0* numOfCloseInterval
    """
    arrayData = array(x)
    arrayData.sort()
    numOfInvalidInterval = 0
    if (len(arrayData)>1):
        for i in xrange(len(arrayData)-1):
            if (arrayData[i+1] - arrayData[i])<0.25:
                numOfInvalidInterval+=1
    return 5.0*numOfInvalidInterval


    
def getBin(x,n):
    """Convert x to n-bit binary string"""
    return x >= 0 and str(bin(x))[2:].zfill(n) or "-" + str(bin(x))[3:].zfill(n)


def t_arrayfit():
    x = [0.236, 0.631, 1.088, 1.515, 2.055, 2.542, 3.173, 3.953,3,510]
    print arrayfit(x)

def t_arraypenalty():
    x = [0.236, 0.631, 1.088, 1.515, 2.055, 2.542, 3.173, 3.253,3,510]
    print arraypenalty(0, x)


if __name__=='__main__':
    use_aopso()
