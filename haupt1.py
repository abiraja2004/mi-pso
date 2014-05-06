#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This program is to verify effectiness of AU-PSO in Low Side-lobe Phase Taper for A Linear Array in 
#    R. Haupt, Antenna design with a mixed integer genetic algorithm,” Antennas Propagation, IEEE Trans., vol. 55, no. 3, pp. 577–582, 2007.

import pylab
import json
from numpy import array, linspace, concatenate, pi, zeros
from arrayfactor import OneDimAF, OneDimSLL
import AOPSO_singleRun as aopso
import RPSO_singleRun as rpso

def optimize_31_elem_array_1():
    #uniform symmetric array
    #phase of centre elem: 0
    #inner 8 elem: real phase
    #next  7 elem: 4-bit phase shifter, precision: pi/8
    def sll(rightphase1):
    # phase[0:8]: real phase
    # phase[8:15]: 4-bit phase shifter
        num=31
        pos=linspace(-15*0.5, 15*0.5, num)
        amp=linspace(1,1,num)
        rightphase = array(rightphase1[:])
        rightphase[8:15] = rightphase[8:15]*(pi/8.0)
        leftphase = rightphase[::-1]
        fullphase = concatenate((leftphase,array([0.]),rightphase))
        AF=zeros(360)
        for i in xrange(360):
            AF[i] = OneDimAF(num, pos, amp, fullphase, i/180.*pi)
        AF = AF - max(AF)
        return OneDimSLL(AF)
    nDim = 15
    numOfParticles = 20
    maxIteration = 300
    minX = concatenate(([0.]*8,[0]*7))
    maxX = concatenate(([2*pi]*8, [16]*7))
    maxV = 0.2*(maxX - minX)
    minV =-1.0*maxV
    step = 3
    numOfRun = 10
    gBestArray = linspace(0,0, 1+maxIteration) 
    f=open('haupt1.txt','w')
    gBest = []
    for i in xrange(numOfRun):
        p = aopso.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, sll, intDim=7, step=step)
        p.run()
        gBest.append(p.gBestArray.tolist())
        gBestArray += p.gBestArray[:]
    json.dump(gBest,f)
    f.close()
    pylab.title('Mixed Array')
    pylab.xlabel('generation')
    pylab.ylabel('max relative sll (dB)')
    pylab.grid(True)
    pylab.plot(range(1+maxIteration), gBestArray/numOfRun)
    pylab.show()

def optimize_31_elem_array_2():
    #uniform symmetric array
    #phase of centre elem: 0
    #all elem: real phase
    def sll(rightphase1):
    # phase[0:8]: real phase
    # phase[8:15]: 4-bit phase shifter
        num=31
        pos=linspace(-15*0.5, 15*0.5, num)
        amp=linspace(1,1,num)
        rightphase = array(rightphase1[:])
        leftphase = rightphase[::-1]
        fullphase = concatenate((leftphase,array([0.]),rightphase))
        AF=zeros(360)
        for i in xrange(360):
            AF[i] = OneDimAF(num, pos, amp, fullphase, i/180.*pi)
        AF = AF - max(AF)
        return OneDimSLL(AF)
    nDim = 15
    numOfParticles = 20
    maxIteration = 300
    minX = array([0.]*15)
    maxX = array([2*pi]*15)
    maxV = 0.2*(maxX - minX)
    minV =-1.0*maxV
    step = 0
    numOfRun = 10
    gBestArray = linspace(0,0, 1+maxIteration) 
    f=open('haupt1_1.txt','w')
    gBest = []
    for i in xrange(numOfRun):
        p = aopso.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, sll, intDim=0, step=step)
        p.run()
        gBest.append(p.gBestArray.tolist())
        gBestArray += p.gBestArray[:]
    json.dump(gBest,f)
    f.close()
    pylab.title('Mixed Array')
    pylab.xlabel('generation')
    pylab.ylabel('max relative sll (dB)')
    pylab.grid(True)
    pylab.plot(range(1+maxIteration), gBestArray/numOfRun)
    pylab.show()



if __name__=='__main__':
    optimize_31_elem_array_2()
