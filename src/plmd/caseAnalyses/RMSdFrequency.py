#!/usr/bin/python
import plmd.plotData as myPlot
import numpy as np

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir ):

    print "Plotting Histogram of RMSd frequencies"
    
    # Get the data
    rmsdValues = []
    with open(caseDir+"/analysis/data/refMinimized.rmsd", "r") as fi:
        next(fi)
        for line in fi:
            rmsdValues.append( float(line.split()[1]) )

    # Create histogram & save in file
    hist, binEdges = np.histogram( rmsdValues , bins=100 , density=True)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    with open(caseDir+"/analysis/data/refMinimizedCleared.rmsd", "w") as fo:
        for i in range(0, len(hist)):
            fo.write( str(bincenters[i])+"\t"+str(hist[i]*100)+"\n" )
        

    print "Max: ",np.max( rmsdValues )
    print "Min: ",np.min( rmsdValues )
    
    # Do the plotting
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "RMSd Frequency", 
        ["Initial Minimization Ref"], 
        [caseDir+"/analysis/data/refMinimizedCleared.rmsd"] , 
        "Frequency (%)",
        xUnit="RMSd ($\AA$)"
    )