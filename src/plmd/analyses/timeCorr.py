#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):
    
    # User info
    print "Creating plot of time correlation for end-to-end distance"    

    # B factor plot
    myPlot.plotData( 
        caseDir , 
        "End-To-End Time Correlation", 
        ["$ \tau $"], 
        ["timeCorr.out"] , 
        "Autocorrelation", 
        skipLines=4)
