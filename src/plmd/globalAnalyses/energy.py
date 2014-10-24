#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( plotTitle, datafiles ):

    print "Creating plot of: "+plotTitle

    myPlot.plotData( 
        "globalAnalysesPlots" , 
        plotTitle, 
        datafiles['caseLabels'], 
        datafiles['filepaths'] , 
        "E (kcal/mol)" ,
        skipLines = 1
    )
    