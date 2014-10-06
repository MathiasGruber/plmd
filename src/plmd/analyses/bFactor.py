#!/usr/bin/python
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # B factor plot
    myPlot.plotData( 
        caseDir , 
        "B Factor for Protein Backbone", 
        ["CA atom", "C atom","N atom", "CA, C and N"], 
        ["bFactor.ca","bFactor.c","bFactor.n","bFactor.all"] , 
        "B Factor", xUnit="Atom ID", types=["bs","gs","rs","-"])