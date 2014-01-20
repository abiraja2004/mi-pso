#!/usr/bin/env python

import AOPSO_singleRun as aopso
import BPSO_singleRun as bpso

def use_aopso():
    



def getBin(x,n):
    #Convert x to n-bit binary string
    return x >= 0 and str(bin(x))[2:].zfill(n) or "-" + str(bin(x))[3:].zfill(n)
