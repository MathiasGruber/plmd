#!/usr/bin/python
import plmd.plotData as myPlot

def runAnalysis( datafiles , resultsDir ):

    print "Creating plot of end-to-end distance"    

    myPlot.plotData( 
        resultsDir+"/plots" , 
        "End to End Residue Distance", 
        datafiles['caseLabels'] , 
        datafiles['filepaths'] , 
        "Distance ($\AA$)",
        scatter = True,
        skipLines = 1
    )
