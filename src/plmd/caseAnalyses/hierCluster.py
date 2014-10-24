#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# A function for separating a dataset based on cluster file
def separateDataSet( dataFile, clusterInfo , outFiles , xColumn = 0 ):
    data = {}
    clusterLen = len(clusterInfo)
    clusters = []
    with open(dataFile, "r") as fi:
        next(fi)
        for line in fi:
            temp = line.split()
            
            # Get the cluster
            if clusterLen >= int(temp[0]):
                
                # Clusters
                clusterNr = clusterInfo[ int(temp[0])-1 ]
                if clusterNr not in clusters:
                    clusters.append(clusterNr)           
                
                # Add this cluster for each component in data file
                for i in range(1, len(temp)):
                    
                    # The identifier for this dataset
                    cluster = "d"+str(i)+"_c"+clusterNr
                    
                    # Check if it already exists
                    if cluster not in data:
                        data[ cluster ] = { "x":[], "y":[] }
                    
                    # Add entry
                    data[ cluster ]["x"].append( temp[ xColumn ] )
                    data[ cluster ]["y"].append( temp[ i ] )

    # Save it in appropriate files
    for cluster, dataset in data.iteritems():
        with open( outFiles+str(cluster) , "w") as fo:
            for i in range( 0, len(dataset["x"]) ):
                fo.write( dataset["x"][i]+"\t"+dataset["y"][i]+"\n" )
      
    return len(clusters)
            

# Function for running the actual analysis
def runAnalysis( caseDir ):

    print "Now Doing hier Cluster"
    
    # Get the cluster vs. frame information
    clusterInfo = []
    with open( caseDir+"/analysis/data/cluster_hier_out" , "r") as fi:
        next(fi)
        for line in fi:
            clusterInfo.append( line.split()[1] )
            
    # RMSd PLOT
    ###########            
            
    # Separate the dataset into separate files.
    numRmsdDataSets = separateDataSet( 
        caseDir+"/analysis/data/backbone.rmsd", 
        clusterInfo ,
        caseDir+"/analysis/data/cluster_hier_rmsd_"
    )
    
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
        "RMSd Cluster Plot", 
        clusterLabels, 
        clusterFiles , 
        "RMSd ($\AA$)",
        scatter = True,
        legendLoc = 4
    )
    
    ## PCA PLOT
    ###########
    
    # Separate the dataset
    numPCAdataSets = separateDataSet( 
        caseDir+"/analysis/data/pca12-ca", 
        clusterInfo ,
        caseDir+"/analysis/data/cluster_hier_pca_",
        xColumn = 1
    )       
    
    # Number of datasets
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
            "PCA Cluster, 1vs"+str(n), 
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