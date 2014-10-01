#!/usr/bin/python

import ConfigParser

config = ConfigParser.RawConfigParser()

# The default config options

config.add_section('simulationParameters')
config.set('simulationParameters', 'forceField', 'leaprc.ff12SB')
config.set('simulationParameters', 'qmCharge', '-2')
config.set('simulationParameters', 'qmShake', '1')
config.set('simulationParameters', 'qmTheory', 'PM6-D')
config.set('simulationParameters', 'ntf', '2')
config.set('simulationParameters', 'ntc', '2')
config.set('simulationParameters', 'timestepSize', '0.002')
config.set('simulationParameters', 'timestepNumber', '1000000')

config.add_section('inputFiles')
config.set('inputFiles', 'ligand', 'filenameLigand.mol2')
config.set('inputFiles', 'ligandCount', '1')
config.set('inputFiles', 'peptide', 'filenamePeptide.pdb')
config.set('inputFiles', 'peptideCount', '1')
config.set('inputFiles', 'cases', '16')

config.add_section('submissionParameters')
config.set('submissionParameters', 'nodeControl', 'nodes=1:ppn=8')
config.set('submissionParameters', 'wallClock', '72:00:00')
config.set('submissionParameters', 'mdRuns', '50')



# Writing our configuration file to 'example.cfg'
with open('defaultExample.cfg', 'wb') as configfile:
    config.write(configfile)