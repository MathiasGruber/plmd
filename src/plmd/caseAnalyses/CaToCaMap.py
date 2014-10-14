#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # User info
    print "Plotting a Ca-to-Ca data map"

    # CA distance map
    myPlot.plotDataMap( 
        caseDir+"/analysis/plots" , 
        "CA Distance Map", 
        caseDir+"/analysis/data/CAdistmat.dat", 
        "Residue ID" , 
        "Residue ID", 
        xColumn=list('SGAGKT'),
        yColumn=list('SGAGKT')
    )
