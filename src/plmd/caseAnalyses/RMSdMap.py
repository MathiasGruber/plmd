#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor ):

    # User info
    print "Plotting a RMSd map"

    # CA distance map
    myPlot.plotDataMap( 
        caseDir+"/analysis/plots" , 
        "RMSd of backbone atoms: C, CA and N", 
        caseDir+"/analysis/data/RMS.ps", 
        "Time (ps)" , 
        "Time (ps)", 
        skipColumn = 1, 
        skipRow = 1, 
        xFactor=timeFactor, 
        yFactor=timeFactor
    )
