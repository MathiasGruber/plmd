#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor ):

     # User info
    print "Creating plot of end-to-end distance"    

    # B factor plot
    myPlot.plotData( 
        caseDir , 
        "End to End Residue Distance", 
        ["End to End"], 
        ["dist_end_to_end.list"] , 
        "Distance ($\AA$)", 
        xFactor = timeFactor )
