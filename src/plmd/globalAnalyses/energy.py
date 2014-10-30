#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( plotTitle, datafiles, resultsDir ):

    print "Creating plot of: "+plotTitle

    myPlot.plotData( 
        resultsDir+"/plots" , 
        plotTitle, 
        datafiles['caseLabels'], 
        datafiles['filepaths'] , 
        "E (kcal/mol)" ,
        skipLines = 1
    )
    