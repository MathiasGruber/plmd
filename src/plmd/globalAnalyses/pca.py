#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs

# Function for running the actual analysis
def runAnalysis( caseDirs ):

    # Determine layout
    rows,columns = 1,1
    dirs = len( caseDirs )
    if dirs > 1:
        columns = 2
        if dirs > 2:
            rows = 2
            
    # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
    # Good URLs for PCA in CPPTRAJ:
    # http://archive.ambermd.org/201404/0243.html

    buffer = "trajin "+caseDirs[0]+"/mergedResult.dcd 1 last 1"
    

    # PCA plotter
    pcaHandler = pcaFuncs.PCA( 
        "globalAnalysesPlots/PCA_analysis.pdf",
        subplotColums = columns,
        subplotRows = rows
    )

    # Go through the case dirs to plot
    for caseDir in caseDirs:

        # Do the plots of energy landscapes & distributions
        pcaHandler.plotPCA( 
            "Case: "+caseDir.split("/")[-1],    # Plot Title
            caseDir+"/analysis/data/" ,         # Data Dir
            "evecs-ca.dat",                     # Eigenvalue file
            "pca12-ca",                         # Eigenvector file
            eigenVectorCount = 2,                # Only plot two
            plotDistibution = False            # Do not plot the distribution
        )
    
    # Save the plot
    pcaHandler.savePlot()
