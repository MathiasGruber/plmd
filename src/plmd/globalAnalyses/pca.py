#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import plmd.generalAnalyses.clusterAnalysis as cluster
import plmd.plotData as myPlot
import os, sys

# Function for running the actual analysis
def runAnalysis( caseDirs , resultsDir ):
    
    # Do a reference for each one
    for refDir in caseDirs:

        # ID of reference case
        refID = refDir.split("/")[-1]
        
        # User info
        print "Doing PCA analysis with "+refDir+" as reference"
        
        # Get the PCA limits of component 1-2 plot
        limit = 10
        with open(refDir+"/analysis/data/pca_limits_1", "r") as fi:
            limit = int(float(fi.read()))
            limit += 0.01
        
        # Go through the case dirs to plot
        for caseDir in caseDirs:
            
            print "Using "+caseDir+" as case"
            
            # ID of case
            caseID = caseDir.split("/")[-1]
            
            # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
            # Good URLs for PCA in CPPTRAJ:
            # http://archive.ambermd.org/201404/0243.html
                        
            # PCA plotter
            pcaHandler = pcaFuncs.PCA( 
                resultsDir+"/plots/pcaComparison/PCA_"+caseID+" on "+refID+".pdf"
            )    
            
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
                "Case: "+caseID+". Ref case: "+refID,   # Plot Title
                caseDir+"/analysis/data/" ,        # Data Dir
                "global_pca",                      # Eigenvector file
                eigenVectorCount = 2,              # Only plot two
                plotDistibution = False,           # Do not plot the distribution
                limits = limit
            )
            
            # Save the plot
            pcaHandler.savePlot()

            # Do both hier and dbscan
            for clusterType in ["dbscan","hier"]:            
                
                # Instantiate the class
                if os.path.isfile(caseDir+"/analysis/data/cluster_"+clusterType+"_out"):   
                    
                    print "Doing the "+clusterType+" cluster equivalent of the PCA plot"
                
                    # Start the cluster handler. Load the file declaring cluster for each frame
                    clusterHandler = cluster.clusterBase( caseDir+"/analysis/data/cluster_"+clusterType+"_out" )
                    
                    # Separate the dataset.
                    # global_pca is the projection file for this case on the ref modes
                    numPCAdataSets = clusterHandler.separateDataSet( 
                        caseDir+"/analysis/data/global_pca",            # Input file
                        caseDir+"/analysis/data/cluster_"+clusterType+"_pca_",   # Output files
                        xColumn = 1
                    ) 
                    
                    # Create lists of labels and files for plotting
                    clusterLabels = []
                    clusterFiles = []
                    offset = 1 if clusterType == "hier" else 0
                    for i in range( 0+offset, numPCAdataSets+offset):
                        clusterLabels.append( "Cluster "+str(i) )
                        clusterFiles.append( caseDir+"/analysis/data/cluster_"+clusterType+"_pca_d2_c"+str(i) )
                    
                    # First one is noise
                    if offset == 0:
                        clusterLabels[0] = "Noise"                 
                    
                    myPlot.plotData( 
                        resultsDir+"/plots/pcaComparison/" , 
                        clusterType+", "+caseID+" on "+refID, 
                        clusterLabels, 
                        clusterFiles , 
                        "PC1",
                        xUnit = "PC1",
                        scatter = True,
                        legendLoc = 4,
                        figWidth = 8,
                        figHeight = 8,
                        tightXlimits = False,
                        legendFrame = 1,
                        legendAlpha = 1,
                        xLimits = [-limit,limit],
                        yLimits = [-limit,limit]
                    )