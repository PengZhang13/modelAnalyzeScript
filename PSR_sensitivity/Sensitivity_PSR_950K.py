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

    print "\n\n\n**************************\n***** Solving case %d *****\n"%(reactionIndex)

    while t < 30.0:
        print t
        t += dt
        network.advance(t)

    results = []
    for spc in targetSpc:
        results.append(cstr.thermo[spc].X[0]*1.0e6)

    end = time.time()
    print 'Execution time is:'
    print end - start

    return [gas.reaction_equation(reactionIndex-1)] + results



def traverseTemperature(TList, P, moleFraction, factor, reactionIndex, targetSpc):

    targetMoleFraction = []
    for T in TList:
        targetMoleFraction.append(runSinglePSR(T, P, moleFraction, factor, reactionIndex, targetSpc))
    with open('output_TList_MoleFraction_30s.csv', 'wb') as f:
        output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
        for i in range(len(TList)):
            output.writerow( [TList[i]]+ targetMoleFraction[i] )


def traverseReactionIndex(T, P, moleFraction, factor, reactionIndexList, targetSpc):

    baseCase = runSinglePSR(T, P, moleFraction, 1.0, 1, targetSpc)

    with open('sensOutput_PSR_' + str(int(T))+ 'K.csv', 'wb') as f:
        output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
        for reactionIndex in reactionIndexList:
            spcData = runSinglePSR(T, P, moleFraction, factor, reactionIndex, targetSpc)
            senDataC3H6 = (spcData[1] - baseCase[1])/baseCase[1]*1.0/(factor-1.0)
            senDataC6H10 = (spcData[2] - baseCase[2])/(8000.0 - baseCase[2])*1.0/(factor-1.0)

            output.writerow( [reactionIndex] + spcData + [senDataC3H6] + [senDataC6H10] )


################################################################################
if __name__ == '__main__':
    start = time.time()

    T = 950.0  # Unit: K
    P = 1.07 *1.0e5    # Unit: bar

    moleFraction = 'C6H10:0.008, HE:0.992'
    factor = 1.25
    reactionIndexList = list(range(1,4527)) # Total number: 1-4526
    targetSpc = ['C3H6', 'C6H10']

    traverseReactionIndex(T, P, moleFraction, factor, reactionIndexList, targetSpc)

#    TList = list(range(800,1075,25))
#    traverseTemperature(TList, P, moleFraction, factor, reactionIndex, targetSpc)
#    print runSinglePSR(T, P, moleFraction, factor, reactionIndex, 'C3H6')

    end = time.time()
    print 'Done! Execution time is:'
    print end - start

