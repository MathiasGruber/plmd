#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    print "Creating plot of potential energy"
    myPlot.plotData( caseDir , "Potential Energy", ["E_pot"], ["summary.EPTOT"] , "E (kcal/mol)" )
    
    print "Creating plot of kinetic energy"    
    myPlot.plotData( caseDir , "Kinetic Energy", ["E_kin"], ["summary.EKTOT"] , "E (kcal/mol)" )
    
    print "Creating plot of total energy"    
    myPlot.plotData( caseDir , "Total Energy", ["E_tot"], ["summary.ETOT"] , "E (kcal/mol)" )