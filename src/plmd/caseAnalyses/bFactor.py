#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # B factor plot
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "B Factor for Protein Backbone", 
        ["CA atom", "C atom","N atom", "CA, C and N"], 
        [caseDir+"/analysis/data/bFactor.ca",caseDir+"/analysis/data/bFactor.c",caseDir+"/analysis/data/bFactor.n",caseDir+"/analysis/data/bFactor.all"] , 
        "B Factor", xUnit="Atom ID", types=["bs","gs","rs","-"])