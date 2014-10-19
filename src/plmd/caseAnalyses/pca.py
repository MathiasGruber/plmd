#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs

# Function for running the actual analysis
def runAnalysis( caseDir, mdTrajectory ):

    # Do the plotting
    pcaHandler = pcaFuncs.PCA( 
        caseDir+"/analysis/plots/PCA_analysis.pdf"
    )

    # Do the plots of energy landscapes & distributions
    pcaHandler.plotPCA( 
        "Case: "+caseDir.split("/")[-1],               # Plot Title
        caseDir+"/analysis/data/" ,                    # Data Dir
        "pca12-ca",                                     # Eigenvector file
        "evecs-ca.dat"                                # Eigenvalue file
    )
    
    # Save all the structures
    pcaHandler.saveStructures( mdTrajectory, caseDir+"/analysis/structures/", caseDir+"/analysis/data/pca12-ca" )
    
    # Save the plot
    pcaHandler.savePlot()
