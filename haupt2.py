#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This program is to optimize the circular polarizated microstrip patch antenna in 
#    R. Haupt, Antenna design with a mixed integer genetic algorithm,” Antennas Propagation, IEEE Trans., vol. 55, no. 3, pp. 577–582, 2007.

import logging
logging.basicConfig(filename='haupt2.log', level=logging.INFO)

import pylab
import os,sys,re
import subprocess
import pdb


from numpy import array, linspace, zeros, pi
import AOPSO_singleRun as aopso


# global variables
fekoFileName='haupt2_parameterized'

def setValue(context, varName, varValue):
    """
    Change variable's value in context
    context -- content of .pre file
    varName -- variable name without leading '#', e.g., 'lam', NOT '#lam'
               if varName does not exist in context, an error occurs
    varValue -- value to be set. type(varValue) should be float
    example: setValue(context, 'lam', 3)
    """
    pattern = '^#'+varName+'=.*\d'
    replace = '#'+varName+'='+str(varValue)
    if (re.search(pattern, context, re.MULTILINE)):
        context = re.sub(pattern, replace, context, flags=re.MULTILINE)
    else:
        sys.exit('.pre file does not contain '+varName)
    return context

def fileExists(filePath):
    return os.path.exists(filePath)

def modifyParameters(preFileName, myDict):
    """Set k,v in myDict to preFileName
    k -- string, variable name
    v -- float, value of k
    """
    if (fileExists(preFileName)):
        context = open(preFileName).read()
        for k,v in myDict.iteritems():
            pattern = '^#'+k+'=.*\d'
            replace = '^#'+k+'='+str(v)
        if (re.search(pattern, context, re.MULTILINE)):
            context = re.sub(pattern, replace, context, flags=re.MULTILINE)
        else:
            sys.exit(preFileName+' does not contain '+k)
        open(preFileName,'w').write(context)


def getSParameter(context, RefZ=50.0):
    """Get S parameter from .out file
    context -- content of .out file
    RefZ -- reference impedance
    context contains only single frequency
    """
    pattern = '^ Impedance in Ohm  .*$'
    Zreal = float(re.search(pattern, context, flags=re.MULTILINE).group().split()[3])
    Zimag = float(re.search(pattern, context, flags=re.MULTILINE).group().split()[4])
    Zl = Zreal+Zimag*1.0j
    Zc = RefZ
    S11 = (Zc-Zl)/(Zc+Zl)
    return abs(S11)

def getE(context):
    """Return  a E(n,4) numpy array
    n -- number of observing points
    E[:,0] -- abs(E_theta)
    E[:,1] -- angle(E_theta)
    E[:,2] -- abs(E_phi)
    E[:,3] -- angle(E_phi)
    """
    pattern = r'angle   direction.*Gain is'  
    Etext = re.search(pattern, context, flags=re.MULTILINE|re.DOTALL).group()
    # Etext = Etext.split('\r\n')[1:-2] #This works on *nix only
    Etext = Etext.split('\n')[1:-2] #This works on Windows only
    n = len(Etext)
    E = zeros((n,4))
    for i in range(n):
        E[i,0] = float(Etext[i].split()[2])
        E[i,1] = float(Etext[i].split()[3])
        E[i,2] = float(Etext[i].split()[4])
        E[i,3] = float(Etext[i].split()[5])
    return E

def haupt_fit(v):
    """Return fitness for CP patch antenna
    v[0] -- Lx=(lambda/4)*(v[0]+1)
    v[1] -- Ly=(lambda/4)*(v[1]+1)
    v[2] -- Fx=0.4*Lx*v[2]
    v[3] -- Fy=0.4*Ly*v[3]
    v[4] -- h=3.15 if v[4]==1, 
            h=1.575 if v[4]==2
    v[5] -- epsi=2.33 if v[5]==1,
            epsi=2.2  if v[5]==2
    """
    logging.info('current v is %s', str(v))
    wavelen = 3e8/10e9*1e3  # Unit: mm
    Lx=wavelen/4*(v[0]+1.0)
    Ly=wavelen/4*(v[1]+1.0)
    Fx=0.4*Lx*v[2]
    Fy=0.4*Ly*v[3]
    if (v[4]==1):
        h=3.15
    elif (v[4]==2):
        h=1.575
    if (v[5]==1):
        myepsi=2.33
    elif (v[5]==2):
        myepsi=2.2
    inName = fekoFileName+'.pre'  #fekoFileName is a global var
    if not(fileExists(inName)):
        sys.exit(inName+' not found')
    context = open(inName).read()
    context = setValue(context, 'Lx', Lx)
    context = setValue(context, 'Ly', Ly)
    context = setValue(context, 'Fx', Fx)
    context = setValue(context, 'Fy', Fy)
    context = setValue(context, 'h', h)
    context = setValue(context, 'myepsi', myepsi)
    logging.info('Lx, Ly, Fx, Fy, h, myepsi = %g, %g, %g, %g, %g, %g',Lx, Ly, 
            Fx, Fy, h, myepsi)
    open(inName,'w').write(context)
    outName = fekoFileName+'.out'
    if fileExists(outName):
        os.remove(outName)
    fekoCmdArgs=[r'C:\Program Files\FEKO\bin\runfeko.exe', inName, '-np', '8']
    runResult = subprocess.check_output(fekoCmdArgs)
    if not('\r\nFinished\r\n' in runResult):
        sys.exit('feko not run correctly')
    if not(fileExists(outName)):
        sys.exit(outName+' not found')
    context = open(outName).read()
    E = getE(context)
    item1 = abs(E[:,0]-E[:,2]) / ((E[:,0]+E[:,2])+1e-7) # to avoid division by zero
    item1 = min(item1)
    item2 = abs(abs(E[:,1]-E[:,3])/180.0*pi - pi/2)
    item2 = min(item2)
    item3 = getSParameter(context)
    return max(array([item1, item2, item3]))

def haupt2():
    logging.info('haupt2 start')
    nDim = 6
    numOfParticles = 12
    maxIteration = 300
    minX = array([0.0, 0.0, 0.0, 0.0, 1, 1])
    maxX = array([1.0, 1.0, 1.0, 1.0, 2, 2])
    maxV = 0.2*(maxX - minX)
    minV = -1.0*(maxX - minX)
    maxV[4], maxV[5] = 1, 1
    minV[4], minV[5] = -1, -1
    step = 0
    numOfRun = 1
    gBestArray = linspace(0, 0, 1+maxIteration)
    gBest = []
    for i in xrange(numOfRun):
        p = aopso.PSOProblem(nDim, numOfParticles, maxIteration, minX, maxX, minV, maxV, haupt_fit, intDim=2, step=step)
        p.run()
        gBestArray += p.gBestArray[:]
    pylab.title('Circular Polarized Patch')
    pylab.xlabel('generation')
    pylab.ylabel('fitness')
    pylab.grid(True)
    pylab.plot(range(1+maxIteration), gBestArray/numOfRun)
    pylab.show()

if __name__=='__main__':
    haupt2()
