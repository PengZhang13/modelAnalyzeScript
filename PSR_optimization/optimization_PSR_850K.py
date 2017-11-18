#! /usr/bin/env python
# encoding:utf-8

"""
Constant volume and solve energy equation.
Brute force sensitivity calculation for a given species.
"""

import sys
import os
import csv
import time
import numpy as np
import cantera as ct
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def runSinglePSR(T, P, moleFraction, factor, reactionIndex, targetSpc):

    start = time.time()

    # create the gas mixture
    gas = ct.Solution('HXD15_Battin_mech.xml')
    gas.TPX = T, P, moleFraction

    # create an upstream reservoir that will supply the reactor. The temperature,
    # pressure, and composition of the upstream reservoir are set to those of the
    # 'gas' object at the time the reservoir is created.
    upstream = ct.Reservoir(gas)

    # Now create the reactor object with the same initial state
    cstr = ct.IdealGasReactor(gas)

    # Set its volume to 80 cm^3. In this problem, the reactor volume is fixed, so
    # the initial volume is the volume at all later times.
    cstr.volume = 80.0*1.0e-6

    # Connect the upstream reservoir to the reactor with a mass flow controller
    # (constant mdot). Set the mass flow rate.
    vdot = 40* 1.0e-6 # m^3/s
    mdot = gas.density * vdot # kg/s
    mfc = ct.MassFlowController(upstream, cstr, mdot=mdot)

    # now create a downstream reservoir to exhaust into.
    downstream = ct.Reservoir(gas)

    # connect the reactor to the downstream reservoir with a valve, and set the
    # coefficient sufficiently large to keep the reactor pressure close to the
    # downstream pressure.
    v = ct.Valve(cstr, downstream, K=1.0e-9)

    # create the network
    network = ct.ReactorNet([cstr])

    # modify the A factor of the given reaction
    gas.set_multiplier(factor,reactionIndex-1)

    # now integrate in time
    t = 0.0
    dt   = 0.1

    while t < 30.0:
        t += dt
        network.advance(t)

    results = []
    for spc in targetSpc:
        results.append(cstr.thermo[spc].X[0]*1.0e6)

    end = time.time()
    print 'Execution time is:'
    print end - start

    return [gas.reaction_equation(reactionIndex-1)] + results



def traverseReactionIndex(T, P, moleFraction, factor, reactionIndexList, targetSpc):

    baseCase = runSinglePSR(T, P, moleFraction, 1.0, 1, targetSpc)

    with open('sensOutput_PSR_' + str(int(T))+ 'K.csv', 'wb') as f:
        output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
        for reactionIndex in reactionIndexList:
            spcData = runSinglePSR(T, P, moleFraction, factor, reactionIndex, targetSpc)
            senDataC3H6 = (spcData[1] - baseCase[1])/baseCase[1]*1.0/(factor-1.0)
            senDataC6H10 = (spcData[2] - baseCase[2])/(8000.0 - baseCase[2])*1.0/(factor-1.0)

            output.writerow( [reactionIndex] + spcData + [senDataC3H6] + [senDataC6H10] )


def traverseFactor(T, P, moleFraction, factorLgList, reactionIndex, targetSpc):

    moleFractionC3H6List = []
    factorList = [10**i for i in factorLgList]
    for i in range(len(factorList)):
        print "\n\n\n**************************\n***** Solving case %d in %d cases*****\n"%(i+1, len(factorList))
        temp = runSinglePSR(T, P, moleFraction, factorList[i], reactionIndex, targetSpc)
        print temp
        moleFractionC3H6List.append(temp[1])

    return factorList, moleFractionC3H6List


def readTargetMoleFraction(targetMoleFractionFile):
    with open(targetMoleFractionFile, 'r') as f:
        stream = f.readlines()

    return float(stream[0])*1.0e6

################################################################################
if __name__ == '__main__':
    start = time.time()

    T = 850.0  # Unit: K
    P = 1.07 *1.0e5    # Unit: bar

    moleFraction = 'C6H10:0.008, HE:0.992'
    reactionIndex = 2280 # C3H5-A + C6H10 = C3H6 + C6H9
    targetSpc = ['C3H6']
    factorLgList = [float(i)/100 for i in range(-50, 10, 10)]
    targetMoleFractionFile = '850.txt'

    [factorList, moleFractionC3H6List] = traverseFactor(T, P, moleFraction, factorLgList, reactionIndex, targetSpc)
    para = UnivariateSpline(moleFractionC3H6List, factorList)

    targetMoleFraction = readTargetMoleFraction(targetMoleFractionFile)
    optimizedFactor = para(targetMoleFraction)
    with open('optimizedFactor.txt', 'wb') as f:
        output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
        output.writerow(['Optimized F:', optimizedFactor])
        output.writerow(['targetMoleFraction', targetMoleFraction])
        for i in range(len(moleFractionC3H6List)):
            output.writerow([factorList[i], moleFractionC3H6List[i]])

    plt.plot(factorList, moleFractionC3H6List, 'r--', optimizedFactor, targetMoleFraction, 'bs')
    plt.xlabel('factor')
    plt.ylabel('C3H6 Mole Fraction / ppm')
    plt.savefig( 'factor_vs_C3H6.png' )
    plt.close()

    end = time.time()
    print 'Done! Execution time is:'
    print end - start

