#! /usr/bin/env python
# encoding:utf-8

"""
Plug flow reactor.
Brute force sensitivity calculation for a given species.
"""

import sys
import os
import csv
import time
import math
import numpy as np
import cantera as ct

def readTemperatureProfile(tempProfile):
    zz = []
    temp = []
    with open(tempProfile, 'r') as f:
        stream = f.readlines()
        for i in range(len(stream)):
            tokens = stream[i].split(',')
            zz.append(float(tokens[0])/100.0)
            temp.append(float(tokens[1][0:-1]))
    return zz, temp


def readReactionIndexFile(reactionIndexFile):

    reactionIndexList = []
    with open(reactionIndexFile, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            reactionIndexList.append(int(row[0]))
    return reactionIndexList


def runSinglePFR(T_0, pressure, composition_0, length, v_0, d, tempProfile, factor, reactionIndex):

    start = time.time()

    # input file containing the reaction mechanism
    reaction_mechanism = 'HXD15_Battin_mech.xml'

    # import the gas model and set the initial conditions
    gas2 = ct.Solution(reaction_mechanism)
    gas2.TPX = T_0, pressure, composition_0
    mass_flow_rate2 = v_0 * gas2.density *(1.01325e5/pressure)*(T_0/298.15)

    # Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
    # of 'n_steps' stirred reactors.
    n_steps = 100
    area = math.pi*d*d/4
    dz = length / n_steps
    r_vol = area * dz

    # create a new reactor
    r2 = ct.IdealGasReactor(gas2)
    r2.volume = r_vol

    # create a reservoir to represent the reactor immediately upstream. Note
    # that the gas object is set already to the state of the upstream reactor
    upstream = ct.Reservoir(gas2, name='upstream')

    # create a reservoir for the reactor to exhaust into. The composition of
    # this reservoir is irrelevant.
    downstream = ct.Reservoir(gas2, name='downstream')

    # The mass flow rate into the reactor will be fixed by using a
    # MassFlowController object.
    m = ct.MassFlowController(upstream, r2, mdot=mass_flow_rate2)

    # We need an outlet to the downstream reservoir. This will determine the
    # pressure in the reactor. The value of K will only affect the transient
    # pressure difference.
    v = ct.PressureController(r2, downstream, master=m, K=1.0e-9)

    sim2 = ct.ReactorNet([r2])

    from scipy.interpolate import UnivariateSpline
    [zz,tem] = readTemperatureProfile(tempProfile)
    para = UnivariateSpline(zz, tem)

    # modify the A factor of the given reaction
    gas2.set_multiplier(factor,reactionIndex-1)

    # iterate through the PFR cells

    # define time, space, and other information vectors
    # z2 = (np.arange(n_steps+1) + 1) * dz
    # t_r2 = np.zeros_like(z2)  # residence time in each reactor
    # u2 = np.zeros_like(z2)
    # t2 = np.zeros_like(z2)

    for n in range(n_steps+1):
        print n
        position = dz*float(n)
        # Set the state of the reservoir to match that of the previous reactor
        gas2.TDY = para(position), r2.thermo.density, r2.thermo.Y

        upstream.syncState()
        sim2.reinitialize()
        # integrate the reactor forward in time until steady state is reached
        sim2.advance_to_steady_state()
        # compute velocity and transform into time
        # u2[n] = mass_flow_rate2 / area / r2.thermo.density
        # t_r2[n] = r2.mass / mass_flow_rate2  # residence time in this reactor
        # t2[n] = np.sum(t_r2)
        # write output data
        # print position, r2.thermo.T, r2.thermo.P, u2[n], t2[n], r2.thermo['C3H6'].X[0]
        # moleFractionC3H6.append([dz*n, u2[n], t2[n], r2.thermo.T, r2.thermo['C3H6'].X[0]])

    end = time.time()
    print 'Done! Execution time is:'
    print end - start

    return [gas2.reaction_equation(reactionIndex-1), r2.thermo['C3H6'].X[0]]



def traverseReactionIndex(T, P, moleFraction, length, vFlowRate, d, tempProfile, factor, reactionIndexList):

    baseCase = runSinglePFR(T, P, moleFraction, length, vFlowRate, d, tempProfile, 1.0, 1)

    with open('sensOutput_PFR_' + tempProfile[0:-4]+ 'K.csv', 'wb') as f:
        output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
        i = 0
        for reactionIndex in reactionIndexList:
            i+=1
            print "\n\n\n**************************\n***** Solving case %d in %d cases*****\n"%(i, len(reactionIndexList))
            spcData = runSinglePFR(T, P, moleFraction, length, vFlowRate, d, tempProfile, factor, reactionIndex)
            senDataC3H6 = (spcData[1] - baseCase[1])/baseCase[1]*1.0/(factor-1.0)

            output.writerow( [reactionIndex] + [spcData[0]] + [spcData[1]] + [senDataC3H6] )




################################################################################
if __name__ == '__main__':

    start = time.time()

    T = 300.0  # inlet temperature [K]
    P = 30.0*133.32  # constant pressure [Pa]
    moleFraction = 'C6H10:0.02, AR:0.98'
    length = 0.6  # *approximate* PFR length [m]
    vFlowRate = 1000.0 * 1.0e-6/60.0 # volumetric flow rate [m3/s]
    d = 6.9e-3 # inner diameter [m]
    tempProfile = '700.csv'
    reactionIndexFile = 'reactionIndexList.csv'

    factor = 1.25
    # reactionIndexList = list(range(1,2)) # Total number: 1-4526
    reactionIndexList = readReactionIndexFile(reactionIndexFile)

    # runSinglePFR(T, P, moleFraction, length, vFlowRate, d, tempProfile)
    traverseReactionIndex(T, P, moleFraction, length, vFlowRate, d, tempProfile, factor, reactionIndexList)

    end = time.time()
    print 'Done! Execution time is:'
    print end - start

