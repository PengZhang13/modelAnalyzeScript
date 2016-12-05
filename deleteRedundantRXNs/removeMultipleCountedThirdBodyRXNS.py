#!/usr/bin/env python
# encoding: utf-8

"""
"""

from rmgpy.chemkin import loadChemkinFile, saveChemkinFile
from rmgpy.rmg.model import ReactionModel
from rmgpy.kinetics.falloff import ThirdBody

# Test if a reaction is a specific form of a third body reaction.
def checkSpecificFormThirdBodyRXN(thirdBodyRXNs, testRNXs):

    M1 = [spc for spc in thirdBodyRXNs.reactants if spc not in testRNXs.reactants]
    M2 = [spc for spc in thirdBodyRXNs.reactants if spc not in testRNXs.products]
    M3 = [spc for spc in thirdBodyRXNs.products if spc not in testRNXs.products]
    M4 = [spc for spc in thirdBodyRXNs.products if spc not in testRNXs.reactants]

    if M1 == [] and M3 == []:
        M5 = [spc for spc in testRNXs.reactants if spc not in thirdBodyRXNs.reactants]
        M6 = [spc for spc in testRNXs.products if spc not in thirdBodyRXNs.products]
        if M5 == M6 or len(M5)*len(M6)==0:
            if checkThirdBodyCoeff(thirdBodyRXNs, [spc for spc in testRNXs.reactants if spc in testRNXs.products][0].molecule[0]):
                return True

    elif M2 == [] and M4 == []:
        M5 = [spc for spc in testRNXs.reactants if spc not in thirdBodyRXNs.products]
        M6 = [spc for spc in testRNXs.products if spc not in thirdBodyRXNs.reactants]
        if M5 == M6 or len(M5)*len(M6)==0:
            if checkThirdBodyCoeff(thirdBodyRXNs, [spc for spc in testRNXs.reactants if spc in testRNXs.products][0].molecule[0]):
                return True

    else:
        return False

# Test
def checkThirdBodyCoeff(thirdBodyRXNs, mol):
    if mol in thirdBodyRXNs.kinetics.efficiencies.keys() and abs(thirdBodyRXNs.kinetics.efficiencies[mol]) <0.001:
        return False
    else:
        return True

def main():

    fuelName = 'test'

    chemFile = 'chem_' + fuelName
    spcFile = 'spc_' + fuelName + '.txt'
    outputChemkinFile = chemFile+'_delete_DupThirdBodyRXNs.inp'
    chemFile += '.inp'

    model = ReactionModel()
    model.species, model.reactions = loadChemkinFile(chemFile, spcFile)
    
    # Delete redundant third body reactions
    thirdBodyRXNs = []
    deleteRedundantThirdBodyRXNs = []
    duplicateEndoRXNs = []

    for rxn in model.reactions:
        if isinstance(rxn.kinetics, ThirdBody):
            thirdBodyRXNs.append(rxn)

    for rxnThirdBody in thirdBodyRXNs:
        for rxnTest in model.reactions:
            if rxnThirdBody.index != rxnTest.index:
                if checkSpecificFormThirdBodyRXN(rxnThirdBody, rxnTest):
                    deleteRedundantThirdBodyRXNs.append(rxnTest)

    # Delete duplicate reactions caused by Intra_R_Add_Endocyclic and Intra_R_Add_Exocyclic
    # Delete the endo one.
    for rxn1 in [rxn for rxn in model.reactions if rxn.duplicate == True and rxn.family == 'Intra_R_Add_Exocyclic']:
        for rxn2 in [rxn for rxn in model.reactions if rxn.duplicate == True and rxn.family == 'Intra_R_Add_Endocyclic']:
            if rxn1.index != rxn2.index and rxn1.isIsomorphic(rxn2, eitherDirection=True):
                duplicateEndoRXNs.append(rxn2)
                rxn1.duplicate = False

    print '--------------------------------------------'
    print 'The following {0} redundant third body reactions will be removed from the original model:'.format(len(deleteRedundantThirdBodyRXNs))
    for rxn in deleteRedundantThirdBodyRXNs:
        print rxn

    print '--------------------------------------------'
    print 'The following {0} duplicated Intra_R_Add_Endocyclic reactions will be removed from the original model:'.format(len(duplicateEndoRXNs))
    for rxn in duplicateEndoRXNs:
        print rxn

    print '--------------------------------------------'
    print 'The original model has {0} rxns'.format(len(model.reactions))
    for rxn in deleteRedundantThirdBodyRXNs + duplicateEndoRXNs:
        model.reactions.remove(rxn)
    print 'The final model has {0} rxns'.format(len(model.reactions))
    
    print 'Save the final model into chemkin file...'
    saveChemkinFile(outputChemkinFile, model.species, model.reactions)

################################################################################
if __name__ == '__main__':
    main()
    
    