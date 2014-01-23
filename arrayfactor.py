#!/usr/bin/env python
"""
Compute array factor of 1d and 2d antenna array.
"""

from numpy import linspace, abs, sum, pi, sin, exp, log10, max, zeros, array,concatenate
import pylab
import ipdb as pdb

#1D array factor
#Unit of AF: dB
def OneDimAF(num, pos, amp, phase, theta):
    # num: number of elements
    # pos: position of elements, with unit of lambda
    # amp: Amplitude of element
    # unit of phase: radian
    # unit of theta: radian
    return 20*log10(abs(sum(amp*exp(1j*(phase/pi+2*pi*(pos-pos[0])*sin(theta))))))

#1D normalized sidelobe level
def OneDimSLL(AF):
    maxGain = max(AF)
#    print 'maxGain is', maxGain
    sll = []
    # find all sidelobe
    for i in range(1,len(AF)-1):
        if AF[i]<maxGain and AF[i-1]<AF[i] and AF[i+1]<AF[i] and (AF[i] not in sll):
            sll.append(AF[i])
    sll.sort()
    return sll[len(sll)-1]
    # return max sidelobe level

def NanbojinTest():
    arrayData = [0.236, 0.631, 1.088, 1.515, 2.055, 2.542, 3.173, 3.953]
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
    print 'SLL of Nanbo Jin''s array is ', OneDimSLL(AF)
#    print 'SLL is ', OneDimSLL(AF)

def NanbojinTest2():
    arrayData = [0.22927538626680652, 0.63628447856160475, 1.0857717937823097, 1.5305810740916364, 2.0523313651720065, 2.5442458800030057, 3.1864080314437637, 3.9516460510870557]
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
    pylab.plot(range(360),AF)
    pylab.show()
    print 'SLL of Nanbo Jin''s array is ', OneDimSLL(AF)
#    print 'SLL is ', OneDimSLL(AF)



def simple_test():
    import pylab
    num = 20
    # distance between element is 1/2*wavelength
    pos = linspace(0,0.5*num,num)
    amp = linspace(1,1,num)
    phase = linspace(0,0,num)
    AF = zeros(360)
    for i in xrange(360):
        AF[i] =  OneDimAF(num, pos, amp, phase, i/180.*pi)
    AF = AF - max(AF)
    pylab.plot(range(360), AF)
    pylab.show()
    print 'SLL of', num, 'elements array is', OneDimSLL(AF), 'dB'

if __name__=='__main__':
    NanbojinTest2()
