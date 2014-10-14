#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    print "Creating plot of potential energy"
    myPlot.plotData( caseDir+"/analysis/plots" , "Potential Energy", ["E_pot"], [caseDir+"/analysis/data/summary.EPTOT"] , "E (kcal/mol)" )
    
    print "Creating plot of kinetic energy"    
    myPlot.plotData( caseDir+"/analysis/plots" , "Kinetic Energy", ["E_kin"], [caseDir+"/analysis/data/summary.EKTOT"] , "E (kcal/mol)" )
    
    print "Creating plot of total energy"    
    myPlot.plotData( caseDir+"/analysis/plots" , "Total Energy", ["E_tot"], [caseDir+"/analysis/data/summary.ETOT"] , "E (kcal/mol)" )