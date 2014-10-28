#!/usr/bin/python
import os
import numpy as np
import plmd.plotData as myPlot

# Function for running the actual analysis
def runAnalysis( caseDirs ):

    # User info
    print "Doing RMSd analysis."

    # Reference Structure
    refStructure = caseDirs[0]+"/pdb-files/finalLEaP_nowat.pdb"

    # For plotting
    dataFiles = []
    dataLabels = []

    # Go through the case dirs to plot
    for caseDir in caseDirs:
        
        # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
        # Good URLs for PCA in CPPTRAJ:
        # http://archive.ambermd.org/201404/0243.html
        
        # Create new submission file
        TEMPLATE = open( caseDir+"/ccptraj_analysis_global_rmsd.ptraj", 'r')
        TEMP = TEMPLATE.read().replace("[REFERENCE]", refStructure  ). \
                               replace("[FOLDER]", caseDir )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseDir+"/ccptraj_analysis_global_rmsd_done.ptraj","w");        
        FILE.write( TEMP );
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+caseDir+"/md-files/peptide_nowat.prmtop -i "+caseDir+"/ccptraj_analysis_global_rmsd_done.ptraj" )

        # Get the data
        rmsdValues = []
        with open(caseDir+"/analysis/data/refGlobalMinimized.rmsd", "r") as fi:
            next(fi)
            for line in fi:
                rmsdValues.append( float(line.split()[1]) )
    
        # Create histogram & save in file
        hist, binEdges = np.histogram( rmsdValues , bins=100 , density=True)
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        with open(caseDir+"/analysis/data/refGlobalMinimizedCleared.rmsd", "w") as fo:
            for i in range(0, len(hist)):
                fo.write( str(bincenters[i])+"\t"+str(hist[i]*100)+"\n" )

        # Save references for plotting
        dataFiles.append( caseDir+"/analysis/data/refGlobalMinimizedCleared.rmsd" )
        dataLabels.append( "Case "+caseDir.split("/")[-1] )

    # Do that plotting    
    myPlot.plotData( 
        "globalAnalysesPlots" , 
        "RMSd Frequency", 
        dataLabels, 
        dataFiles , 
        "Frequency (%)",
        xUnit="RMSd ($\AA$)",
        skipLines = 1
    )
    