#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import plmd.plotData as myPlot

# Function for running the actual analysis
def runAnalysis( caseDir, mdTrajectory ):

    # Create eigenvalue file
    with open( caseDir+"/analysis/data/evecs-cleaned.dat", "w" ) as fo:
        with open(caseDir+"/analysis/data/evecs-ca.dat","r") as f:
            lookForVec = 1
            for line in f:
                temp = line.split()
                if len(temp) == 2 and temp[0] == str(lookForVec):
                    if lookForVec < 15:
                        fo.write( str(lookForVec)+"\t"+temp[1]+"\n" )
                        lookForVec += 1

    # Do the plotting
    pcaHandler = pcaFuncs.PCA( 
        caseDir+"/analysis/plots/PCA_analysis.pdf"
    )

    # Do the plots of energy landscapes & distributions
    pcaHandler.plotPCA( 
        "Case: "+caseDir.split("/")[-1],              # Plot Title
        caseDir+"/analysis/data/" ,                   # Data Dir
        "pca12-ca",                                   # Eigenvector file
        "evecs-cleaned.dat"                                # Eigenvalue file
    )
    
    # Save all the structures
    pcaHandler.saveStructures( mdTrajectory, caseDir+"/analysis/structures/", caseDir+"/analysis/data/pca12-ca" )
    
    # Save the plot
    pcaHandler.savePlot()
    
    # Save an eigenvalue vs eigenvector plot
    print "Creating plot of eigenvalues"
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "Eigenvalue Plot", 
        ["Eigenvalue"], 
        [caseDir+"/analysis/data/evecs-cleaned.dat"] , 
        "AU",
        scatter = True
    )
    
