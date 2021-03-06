

#!/usr/bin/python
import ConfigParser

# Get a default config file
def getDefaultConfig():
    
    # The default config options
    config = ConfigParser.RawConfigParser()
    
    config.add_section('plmd_settings') 
    config.set('plmd_settings', 'name', 'SimulationIdentifier')
    
    config.add_section('inputFiles')
    config.set('inputFiles', 'ligand', 'filenameLigand.mol2')
    config.set('inputFiles', 'ligandCount', '1')
    config.set('inputFiles', 'pLigand', 'false')
    config.set('inputFiles', 'peptide', 'filenamePeptide.pdb')
    config.set('inputFiles', 'peptideCount', '1')
    config.set('inputFiles', 'pPeptide', 'false')
    config.set('inputFiles', 'cases', '4')
    config.set('inputFiles', 'quiet', 'false')
    config.set('inputFiles', 'noTranslate', 'false')
    config.set('inputFiles', 'waterboxsize', '8.0')
    
    config.add_section('mmpbsaSettings')
    config.set('mmpbsaSettings', 'interval', '1000')
    
    config.add_section('simulationParameters')
    config.set('simulationParameters', 'forceField', 'leaprc.ff14SB')
    config.set('simulationParameters', 'polarizeable', 'false')
    config.set('simulationParameters', 'qmEnable', 'true')
    config.set('simulationParameters', 'qmCharge', '-2')
    config.set('simulationParameters', 'qmShake', '1')
    config.set('simulationParameters', 'qmTheory', 'PM6-D')
    config.set('simulationParameters', 'qmRegionOverwrite', 'false')
    config.set('simulationParameters', 'ntf', '2')
    config.set('simulationParameters', 'ntc', '2')
    config.set('simulationParameters', 'ntt', '3')
    config.set('simulationParameters', 'ntb', '2')
    config.set('simulationParameters', 'ntp', '1')
    config.set('simulationParameters', 'gamma_ln', '2.0')
    config.set('simulationParameters', 'timestepPerFrame', '1000')
    config.set('simulationParameters', 'timestepSize', '0.002')
    config.set('simulationParameters', 'timestepNumber', '1000000')

    config.add_section('submissionParameters')
    config.set('submissionParameters', 'nodeControl', 'nodes=2:ppn=8')
    config.set('submissionParameters', 'wallClock', '72:00:00')
    config.set('submissionParameters', 'mdRuns', '50')

    config.add_section('analysisParameters')
    config.set('analysisParameters', 'noMerge', 'false')
    config.set('analysisParameters', 'noStrip', 'false')
    config.set('analysisParameters', 'noBlock', 'false')
    config.set('analysisParameters', 'noEnergy', 'false')
    config.set('analysisParameters', 'noReweight', 'false')
    config.set('analysisParameters', 'noFTP', 'false')
    config.set('analysisParameters', 'dbscanEps', '1.2')
    config.set('analysisParameters', 'dbscanMinPoints', '20')
    config.set('analysisParameters', 'clusterFrames', '10000')
    
    config.add_section('gpuRun')
    config.set('gpuRun', 'enable', 'false')
    config.set('gpuRun', 'gpus', '1')
    
    config.add_section('amdRun')
    config.set('amdRun', 'enable', 'false')
    config.set('amdRun', 'iamd', '2')
    config.set('amdRun', 'energyPerResidue', '4')
    config.set('amdRun', 'energyPerAtom', '0.16')
    config.set('amdRun', 'alphaDfactor', '0.2')
    
    config.add_section('ftpConfiguration')
    config.set('ftpConfiguration', 'server', 'ftp.nanomathias.com')
    config.set('ftpConfiguration', 'ftpUser', 'plmd')
    
    return config