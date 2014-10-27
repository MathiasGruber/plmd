#!/usr/bin/python
import plmd.plotData as myPlot
import plmd.generalAnalyses.clusterAnalysis as cluster

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # User information
    print "Now Doing hier Cluster"
    
    # Instantiate the class
    handler = cluster.clusterBase( caseDir+"/analysis/data/cluster_hier_out" )    
      
    # RMSd PLOT
    numRmsdDataSets = handler.separateDataSet( 
        caseDir+"/analysis/data/backbone.rmsd",
        caseDir+"/analysis/data/cluster_hier_rmsd_"
    )
    
    # User info
    print "Number of datasets for RMSd: "+str(numRmsdDataSets) 
    
    # Create lists of labels and files for plotting
    clusterLabels = []
    clusterFiles = []
    for i in range( 1, numRmsdDataSets+1):
        clusterLabels.append( "Cluster "+str(i) )
        clusterFiles.append( caseDir+"/analysis/data/cluster_hier_rmsd_d1_c"+str(i) )
    
    # Do the plottin
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "RMSd Hier Cluster Plot", 
        clusterLabels, 
        clusterFiles , 
        "RMSd ($\AA$)",
        scatter = True,
        legendLoc = 4
    )
    
    ## PCA PLOT
    ###########
    
    # Separate the dataset
    numPCAdataSets = handler.separateDataSet( 
        caseDir+"/analysis/data/pca12-ca", 
        caseDir+"/analysis/data/cluster_hier_pca_",
        xColumn = 1
    )       
    
    # User info
    print "Number of datasets for PCA: "+str(numPCAdataSets)    
    
    # Go through each PCA component
    for n in range(2,4):    
    
        # Create lists of labels and files for plotting
        clusterLabels = []
        clusterFiles = []
        for i in range( 1, numPCAdataSets+1):
            clusterLabels.append( "Cluster "+str(i) )
            clusterFiles.append( caseDir+"/analysis/data/cluster_hier_pca_d"+str(n)+"_c"+str(i) )
        
        print "Calling plot"
        print "labels", clusterLabels  
        print "files", clusterFiles  
        
        myPlot.plotData( 
            caseDir+"/analysis/plots" , 
            "PCA Hier Cluster, 1vs"+str(n), 
            clusterLabels, 
            clusterFiles , 
            "PC"+str(n),
            xUnit = "PC1",
            scatter = True,
            legendLoc = 4,
            figWidth = 8,
            figHeight = 8,
            tightXlimits = False,
            legendFrame = 1,
            legendAlpha = 1
        )