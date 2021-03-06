#!/usr/bin/python
import os, sys
import numpy as np
import plmd.plotData as myPlot

# Function for running the actual analysis
def runAnalysis( caseDirs , resultsDir , noReweight = False ):

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
        TEMPLATE = open( caseDir+"/ccptraj_analysis_rmsd.ptraj", 'r')
        TEMP = TEMPLATE.read().replace("[REFERENCE]", refStructure  ). \
                               replace("[FOLDER]", caseDir )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseDir+"/ccptraj_analysis_rmsd_done.ptraj","w");        
        FILE.write( TEMP );
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+caseDir+"/md-files/peptide_nowat.prmtop -i "+caseDir+"/ccptraj_analysis_rmsd_done.ptraj" )

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
        
        # Check if we should do a reweighted version
        if noReweight == False:
            if os.path.isfile( caseDir+"/md-logs/weights.dat" ):
                
                # User info
                print "aMD weights found. Now attempting reweighting"
                
                # Prepare input file
                numLines = 0
                with open(caseDir+"/analysis/data/refGlobalMinimized.rmsd", "r") as fi:
                    with open(caseDir+"/analysis/data/refGlobalMinimized_singleColumn.rmsd", "w") as fo:
                        next(fi)
                        for line in fi:
                            numLines += 1
                            fo.write( line.split()[1]+"\n" )
                
                # Run the reweighting procedure
                os.system("python $PLMDHOME/src/PyReweighting/PyReweighting-1D.py \
                            -input "+caseDir+"/analysis/data/refGlobalMinimized_singleColumn.rmsd \
                            -name "+caseDir+"/analysis/data/refGlobalMinimized_reweighted \
                            -Xdim 0 6 \
                            -disc 0.06 \
                            -cutoff 5 \
                            -Emax 200 \
                            -job amdweight_CE \
                            -weight "+caseDir+"/md-logs/weights.dat | tee -a reweight_variable.log")
                            
                # Save references for plotting            
                dataFiles.append( caseDir+"/analysis/data/refGlobalMinimized_reweighted-hist-c2.dat" )
                dataLabels.append( "Case "+caseDir.split("/")[-1]+" Reweigthed 2" )
            
    # Do that plotting    
    myPlot.plotData( 
        resultsDir+"/plots" , 
        "RMSd Frequency", 
        dataLabels, 
        dataFiles , 
        "Frequency (%)",
        xUnit="RMSd ($\AA$)",
        skipLines = 1
    )
    