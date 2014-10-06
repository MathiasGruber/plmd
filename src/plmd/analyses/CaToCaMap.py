#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # Plotting a Ca-to-Ca data map

    # CA distance map
    myPlot.plotDataMap( 
        caseDir, 
        "CA Distance Map", 
        "CAdistmat.dat", 
        "Residue ID" , 
        "Residue ID", 
        xColumn=list('SGAGKT'),
        yColumn=list('SGAGKT')
    )
