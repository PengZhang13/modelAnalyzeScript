#!/usr/bin/env python
# encoding: utf-8

"""
This script uses Chemkin to run simulations of a rapid compression machine 
using an RMG-generated model. 
"""

import os
import csv
import math
import numpy as np
from numpy import polyval
import matplotlib.pyplot as plt
#from numpy.polynomial.polynomial import polyval

def readTargetTime(targetTimeFile):
    targetTimeList = []
    targetMoleListC3H6 = []
    with open(targetTimeFile, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            targetTimeList.append(float(row[0])/1000)
            targetMoleListC3H6.append(float(row[1]))
    return targetTimeList, targetMoleListC3H6


def extractTargetData(outputDir, targetTimeFile):
    
    targetTimeList, targetMoleListC3H6 = readTargetTime(targetTimeFile)
    targetData = []
    
    for root, dirs, files in os.walk(outputDir):
        for filename in files:
            filename = os.path.join(outputDir, filename)
            targetData.append(processSingleFile(filename, targetTimeList))
    
    # print targetData[0][1][1]
    
    polyfitCoeffC3H6 = []
    res = []
    for i in range(len(targetData[0][1])):
        fList = []
        C3H6List = []
        HXD15List = []
        
        fittingDir = 'fitting/'
        if not os.path.exists(fittingDir):
            os.makedirs(fittingDir)
        
        fname = fittingDir + 'case6_t_'+str(targetData[0][1][i]*1000)+'ms'
        
        with open(fname + '.csv', 'wb') as f:
            output = csv.writer(f, delimiter = ',', quoting = csv.QUOTE_ALL)
            for j in range(len(targetData)):
                output.writerow([10**float(targetData[j][0]), targetData[j][2][i]*1e6, targetData[j][3][i]*1e6])

                fList.append(10**float(targetData[j][0]))
                C3H6List.append(targetData[j][2][i]*1e6)
                HXD15List.append(targetData[j][3][i]*1e6)

        aaa = np.polyfit( np.array(fList), np.array(C3H6List), 3, full = True)
        polyfitCoeffC3H6.append( aaa[0] )
        res.append(aaa[1][0])

        x_f = [10**(float(i)/100) for i in range(-100,101)]
        y_f = polyval( aaa[0], x_f )
        plt.plot(x_f, y_f, 'r--', np.array(fList), np.array(C3H6List), 'bs')
        plt.xlabel('factor')
        plt.ylabel('C3H6 Mole Fraction / ppm')
        plt.savefig( fname + 'FittingCurve.png' )
        plt.close()

    optF = optimizeFactor(targetTimeList, targetMoleListC3H6, polyfitCoeffC3H6)

    filename_factor = []
    for root, dirs, files in os.walk(outputDir):
        for filename in files:
            filename = os.path.join(outputDir, filename)
            token = filename.split('_F_')
            filename_factor.append(token[1][0:-4])
    bbb = [abs(float(i) - optF) for i in filename_factor]
    optFactorFile = 'output/output_F_' + filename_factor[bbb.index(min(bbb))] + '.csv'
    optTime = []
    optC3H6 = []
    with open(optFactorFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            optTime.append(row[0])
            optC3H6.append(row[1])
    optTime = [float(i)*1000-40 for i in optTime]
    optC3H6 = [float(i)*1e6 for i in optC3H6]
    targetTimeList = [i*1000 for i in targetTimeList]
    plt.plot(optTime, optC3H6, 'r--', targetTimeList, targetMoleListC3H6, 'bs')
    plt.xlabel('Time / ms')
    plt.ylabel('C3H6 Mole Fraction / ppm')
    plt.savefig( 'optimized_C3H6_profile.png' )
    plt.close()


def optimizeFactor(timeList, targetList, polyfitCoeffList):

    fList = []
    resList = []
    for f in range(-100000,100000):
        res = sum([(polyval( polyfitCoeffList[i], 10**(float(f)/100000) ) - targetList[i])**2 for i in range(len(timeList))])
        fList.append(10**(float(f)/100000))
        resList.append(res)
    minIndex = resList.index(min(resList))

    print minIndex, fList[minIndex], resList[minIndex]
    with open('optimized_parameter.txt', 'wb') as stream:
        output = csv.writer(stream, delimiter = ',', quoting = csv.QUOTE_ALL)
        output.writerow(['Optimized F:', fList[minIndex]])
        output.writerow(['Residual:', resList[minIndex]])

    return math.log10(fList[minIndex])

def processSingleFile(filename, targetTimeList):
    token = filename.split('_F_')
    factor = token[1][0:-4]
    
    time = []
    C3H6 = []
    HXD15 = []

    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            time.append(float(row[0])-0.04)
            C3H6.append(float(row[1]))
            HXD15.append(float(row[2]))

    targetTimeTime = []
    targetTimeC3H6 = []
    targetTimeHXD15 = []

    for targetTime in targetTimeList:
        for i in range(len(time)):
            if time[i]<targetTime and time[i+1]>targetTime:
                ratio = (targetTime - time[i])/(time[i+1] - time[i])
                targetTimeC3H6.append(interpolate(ratio, C3H6[i], C3H6[i+1]))
                targetTimeHXD15.append(interpolate(ratio, HXD15[i], HXD15[i+1]))
                targetTimeTime.append(targetTime)
                break

    return [factor, targetTimeTime, targetTimeC3H6, targetTimeHXD15]


def interpolate(ratio, a, b):
    return (b-a)*ratio+a


################################################################################
if __name__ == '__main__':

    currentDir = os.path.dirname(__file__)
    outputDir = os.path.join(currentDir, 'output')
    targetTimeFile = os.path.join(currentDir, 'case6_C3H6.csv')

    extractTargetData(outputDir, targetTimeFile)


