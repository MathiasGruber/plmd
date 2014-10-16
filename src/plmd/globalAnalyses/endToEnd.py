#!/usr/bin/python
import plmd.plotData as myPlot

def runAnalysis( datafiles ):

    print "Creating plot of end-to-end distance"    

    myPlot.plotData( 
        "globalAnalysesPlots" , 
        "End to End Residue Distance", 
        datafiles['caseLabels'] , 
        datafiles['filepaths'] , 
        "Distance ($\AA$)"
    )
