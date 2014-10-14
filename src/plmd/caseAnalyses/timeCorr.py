#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor ):
    
    # User info
    print "Creating plot of time correlation for end-to-end distance"    

    # B factor plot
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "End-To-End Time Correlation", 
        ["$ \tau $"], 
        [caseDir+"/analysis/data/timeCorr.out"] , 
        "Autocorrelation", 
        skipLines=4,
        xFactor = timeFactor )
