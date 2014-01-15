# Test functions

from numpy import cos,pi, sqrt

def Sphere(x):
    result = 0.0
    for i in xrange(len(x)):
        result += x[i]**2 
    return result

def Rosenbrock(x):
    result = 0.0
    for i in xrange(len(x)-1):
        result += (1.0-x[i])**2 + 100.0*((x[i+1]-x[i]**2)**2)
    return result

def Rastrigrin(x):
    result = 0.0
    for i in xrange(len(x)):
        result += x[i]**2 - 10.0*cos(2*pi*x[i])+10
    return result

def Griewank(x):
    result = 0.0
    part1 = 0.0
    part2 = 1.0
    for i in xrange(len(x)):
        part1 += x[i]*x[i]
        part2 *= cos(x[i]/sqrt(i+1))
    result = 1.0/4000*part1-part2+1
    return result
