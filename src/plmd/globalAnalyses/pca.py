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

    # Do the plotting
    pcaHandler = pcaFuncs.PCA( 
        "globalAnalysesPlots/PCA_analysis.pdf",
        subplotColums = columns,
        subplotRows = rows
    )

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
