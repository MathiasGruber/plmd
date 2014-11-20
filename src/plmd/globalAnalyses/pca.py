#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import os, sys

# Function for running the actual analysis
def runAnalysis( caseDirs , resultsDir ):
    
    # User info
    print "Doing PCA analysis."
    
    # Determine layout
    rows,columns = 1,1
    dirs = len( caseDirs )
    if dirs > 1:
        columns = 2
        if dirs > 2:
            rows = 2
            
    
    # Do a reference for each one
    for refDir in caseDirs:

        # ID of reference case
        refID = refDir.split("/")[-1]
        
        # Get the PCA limits of component 1-2 plot
        limit = 10
        with open(refDir+"/analysis/data/pca_limits_1", "r") as fi:
            limit = int(fi.read())
        
         # PCA plotter
        pcaHandler = pcaFuncs.PCA( 
            resultsDir+"/plots/PCA_analysis_Components"+refID+".pdf",
            subplotColums = columns,
            subplotRows = rows
        )

        # Go through the case dirs to plot
        for caseDir in caseDirs:
            
            # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
            # Good URLs for PCA in CPPTRAJ:
            # http://archive.ambermd.org/201404/0243.html
            
            # Create new submission file
            TEMPLATE = open( caseDir+"/ccptraj_analysis_pca.ptraj", 'r')
            TEMP = TEMPLATE.read().replace("[PCAREFERENCE]", refDir  )
            TEMPLATE.close()
                                  
            # Write the submission file
            FILE = open(caseDir+"/ccptraj_analysis_pca.ptraj","w");        
            FILE.write( TEMP );
            FILE.close();
            
            # Run the cpptraj utility
            os.system( "$AMBERHOME/bin/cpptraj -p "+caseDir+"/md-files/peptide_nowat.prmtop -i "+caseDir+"/ccptraj_analysis_pca.ptraj" )
        
    
            # Do the plots of energy landscapes & distributions
            pcaHandler.plotPCA( 
                "Case: "+caseDir.split("/")[-1]+". Ref case: "+refID,   # Plot Title
                caseDir+"/analysis/data/" ,        # Data Dir
                "global_pca",                      # Eigenvector file
                eigenVectorCount = 2,              # Only plot two
                plotDistibution = False,           # Do not plot the distribution
                limits = limit
            )
            
        # Save the plot
        pcaHandler.savePlot()