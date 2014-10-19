#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import os

# Function for running the actual analysis
def runAnalysis( caseDirs ):

    # Determine layout
    rows,columns = 1,1
    dirs = len( caseDirs )
    if dirs > 1:
        columns = 2
        if dirs > 2:
            rows = 2
            
    # PCA plotter
    pcaHandler = pcaFuncs.PCA( 
        "globalAnalysesPlots/PCA_analysis.pdf",
        subplotColums = columns,
        subplotRows = rows
    )

    # Go through the case dirs to plot
    for caseDir in caseDirs:
        
        # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
        # Good URLs for PCA in CPPTRAJ:
        # http://archive.ambermd.org/201404/0243.html
        
        # Create new submission file
        TEMPLATE = open( caseDir+"/ccptraj_analysis_global.ptraj", 'r')
        TEMP = TEMPLATE.read().replace("[PCAREFERENCE]", caseDirs[0]  )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseDir+"/ccptraj_analysis_global_done.ptraj","w");        
        FILE.write( TEMP );
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+caseDir+"/md-files/peptide_nowat.prmtop -i "+caseDir+"/ccptraj_analysis_global_done.ptraj" )
    

        # Do the plots of energy landscapes & distributions
        pcaHandler.plotPCA( 
            "Case: "+caseDir.split("/")[-1],   # Plot Title
            caseDir+"/analysis/data/" ,        # Data Dir
            "global_pca",                      # Eigenvector file
            eigenVectorCount = 2,              # Only plot two
            plotDistibution = False            # Do not plot the distribution
        )
    
    # Save the plot
    pcaHandler.savePlot()
