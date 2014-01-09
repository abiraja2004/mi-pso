"""
Compute array factor of 1d and 2d antenna array.
"""

from numpy import linspace, abs, sum, pi, sin, exp, log10, max, zeros

def ArrayFactor(num, pos, amp, phase, theta):
    # num: number of elements
    # pos: position of elements, with unit of lambda
    # amp: Amplitude of element
    # unit of phase: radian
    # unit of theta: radian
    return 20*log10(abs(sum(exp(1j*(phase/pi+2*pi*(pos-pos[0])*sin(theta))))))

def SLL(AF):
    maxGain = max(AF)
#    print 'maxGain is', maxGain
    sll = []
    # find all sidelobe
    for i in range(1,len(AF)-1):
        if AF[i]<maxGain and AF[i-1]<AF[i] and AF[i+1]<AF[i] and (AF[i] not in sll):
            sll.append(AF[i])
    sll.sort()
    return sll[len(sll)-2]
    # return max sidelobe level
    return max(a1)

if __name__=='__main__':
    import pylab
    num = 20
    # distance between element is 1/2*wavelength
    pos = linspace(0,0.5*num,num)
    amp = linspace(1,1,num)
    phase = linspace(0,0,num)
    AF = zeros(360)
    for i in xrange(360):
        AF[i] =  ArrayFactor(num, pos, amp, phase, i/180.*pi)
    AF = AF - max(AF)
    pylab.plot(range(360), AF)
    pylab.show()
    print 'SLL of', num, 'elements array is', SLL(AF), 'dB'
