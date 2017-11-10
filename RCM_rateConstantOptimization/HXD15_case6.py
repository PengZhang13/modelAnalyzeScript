#!/usr/bin/env python
# encoding: utf-8

"""
This script uses Chemkin to run simulations of a rapid compression machine 
using an RMG-generated model. 
"""

import os
import csv
import numpy

from pychemkin.chemkin import ChemkinJob
from pychemkin.chemkin import getVariable

def runHXD15Job(chemFileName, tempDirName, outputFile):
    import time
    start = time.time()
    job = ChemkinJob(name = 'HXD15_pychemkin', chemFile = chemFileName, tempDir = tempDirName)

# Solver time step profile
#     solTimeStepFile=os.path.realpath(os.path.join(currentDir, 'case1_t.csv'))

    job.preprocess()

    T = 304.15      # Temperature (K)
    P = 1.52944     # Pressure (bar)

    HXD15 = 0.003004
    ar = 1.0 - HXD15
    endtime = 0.055

    input = job.writeInputHomogeneousBatch(
    
        problemType = 'constrainVandSolveE', # Problem type(constrainVandSolveE,
                                           # constrainPandSolveE, constrainPandT, constrainVandT)

        reactants=[('C6H10',HXD15),
                 ('AR',ar),],     # Reactant (mole fraction)

        temperature = T, # Temperature(K)
        pressure = P ,   # Pressure (bar)
        endTime = endtime ,   # End Time (sec)

        #Continuations = True,             # Continuations
        #typeContinuation = 'NEWRUN',      # Type of continuation NEWRUN or CNTN
        #Tlist=[1000.,1200.],              # Temperature (K) list of continuations
        #Plist=[1.,2.],                   # Pressure (atm) list of continuations

        variableVolume = True,            # Variable volume true / false
        variableVolumeProfile ='case6_V.csv',  # File name path for vpro (sec,cm3)

        solverMaxTimeStep = True,
        maxTimeStep = 1.0E-4,
        # solverTimeStepProfile = solTimeStepFile # Solver time profile (sec)
        )
        
    job.run(input, model='CKReactorGenericClosed', pro=True)
    print 'Finish running!\n'
    job.postprocess(sens=False, rop=False, all=True, transpose=False)
    print 'Finish postprocessing!\n'

    print 'CHEMKIN JOBS DONE!'
    
    extractCKCSVFile(job.ckcsvFile, outputFile)

    end = time.time()
    print 'Execution time is:'
    print end - start


def getMoleFraction(ckcsvFile, species=[]):
# The getMoleFraction function from pychemkin is unsuitable in this case. So I decided to make a new one.
    spcData = []
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if label.startswith('Mole_fraction') and len(tokens) == 3:
                for spec in species:
                    if tokens[2] == spec:
                        spcData.append([float(r) for r in row[2:]])
    return spcData

def extractCKCSVFile(ckcsvFile, outputFile):

    [timeData, pressureData] = getVariable(ckcsvFile, 'pressure')
    [timeData, temperatureData] = getVariable(ckcsvFile, 'temperature')
    [timeData, volumeData] = getVariable(ckcsvFile, 'volume')

    spcData = getMoleFraction(ckcsvFile, ['C3H6','C6H10'])
    # In the HXD15 case, only Mole fractions for these species are needed.

    with open(outputFile, 'wb') as f:
        output = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)
        for i in range(len(timeData[0])):
            output.writerow([timeData[0][i], spcData[0][i], spcData[1][i]])
            # output.writerow([timeData[0][i], pressureData[0][i], temperatureData[0][i], volumeData[0][i], spcData[0][i], spcData[1][i]])

    print 'Data has been extracted from the CKCSV file!'


def modifySingleInp(inputFile, modifyInpFile, factor):
    
    with open(inputFile, 'r') as f:
        stream = f.readlines()
        tokens = stream[8322].split('\t')
        # print len(tokens)
        # for i in range(len(tokens)):
        #     print tokens[i]
        tokens[1] = str(float(tokens[1])*factor)
        stream[8322] = '\t'.join(tokens)
        print '\nThe modified line is:'
        print stream[8322]

    with open(modifyInpFile, 'w') as f:
        f.writelines(stream)


################################################################################
if __name__ == '__main__':

    import time
    start = time.time()
    
    currentDir = os.path.dirname(__file__)
    
    inputFile = os.path.join(currentDir, 'HXD15_Battin_mech_therm_base.inp')
    modifyInpFile = os.path.join(currentDir, inputFile[0:-4]+'_mod.inp')
    tempDir = os.path.join(currentDir, 'temp')
    outputDir = os.path.join(currentDir, 'output')
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    factor = [float(i)/100 for i in range(-56, 102, 2)]

    i = 0
    for f in factor:
        print '\n\n\n*********************************\n*******CALCULATION STARTED*******\n*******CASE %d IN %d CASES*******'%(i+1, len(factor))
        i+=1
        outputFile = os.path.join(outputDir, 'output_F_')+str(f)+'.csv'
        
        modifySingleInp(inputFile, modifyInpFile, 10**f)
        runHXD15Job(modifyInpFile, tempDir, outputFile)

    end = time.time()
    print 'Done! Execution time is:'
    print end - start


