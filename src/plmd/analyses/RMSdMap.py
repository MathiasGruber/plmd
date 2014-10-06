#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor ):

    # User info
    print "Plotting a RMSd map"

    # CA distance map
    myPlot.plotDataMap( 
        caseDir, 
        "RMSd of backbone atoms: C, CA and N", 
        "RMS.ps", 
        "Time (ps)" , 
        "Time (ps)", 
        skipColumn = 1, 
        skipRow = 1, 
        xFactor=timeFactor, 
        yFactor=timeFactor
    )
