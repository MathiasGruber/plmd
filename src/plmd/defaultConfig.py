#!/usr/bin/python
import ConfigParser

# Get a default config file
def getDefaultConfig():
    
    # The default config options
    config = ConfigParser.RawConfigParser()
    
    config.add_section('inputFiles')
    config.set('inputFiles', 'ligand', 'filenameLigand.mol2')
    config.set('inputFiles', 'ligandCount', '1')
    config.set('inputFiles', 'pLigand', 'false')
    config.set('inputFiles', 'peptide', 'filenamePeptide.pdb')
    config.set('inputFiles', 'peptideCount', '1')
    config.set('inputFiles', 'pPeptide', 'false')
    config.set('inputFiles', 'cases', '10')
    config.set('inputFiles', 'quiet', 'false')
    
    config.add_section('simulationParameters')
    config.set('simulationParameters', 'forceField', 'leaprc.ff12SB')
    config.set('simulationParameters', 'qmCharge', '-2')
    config.set('simulationParameters', 'qmShake', '1')
    config.set('simulationParameters', 'qmTheory', 'PM6-D')
    config.set('simulationParameters', 'ntf', '2')
    config.set('simulationParameters', 'ntc', '2')
    config.set('simulationParameters', 'timestepSize', '0.002')
    config.set('simulationParameters', 'timestepNumber', '1000000')

    config.add_section('submissionParameters')
    config.set('submissionParameters', 'nodeControl', 'nodes=1:ppn=8')
    config.set('submissionParameters', 'wallClock', '72:00:00')
    config.set('submissionParameters', 'mdRuns', '50')

    config.add_section('analysisParameters')
    config.set('analysisParameters', 'email', 'nano.mathias@gmail.com')
    config.set('analysisParameters', 'noMerge', 'false')
    config.set('analysisParameters', 'noStrip', 'false')
    
    return config