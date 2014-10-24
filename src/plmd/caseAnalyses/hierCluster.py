#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# A function for separating a dataset based on cluster file
def separateDataSet( dataFile, clusterInfo , outFiles ):
    data = {}
    clusterLen = len(clusterInfo)
    with open(dataFile, "r") as fi:
        next(fi)
        for line in fi:
            temp = line.split()
            
            # Get the cluster
            if clusterLen >= int(temp[0]):
                cluster = clusterInfo[ int(temp[0])-1 ]
                
                # Check if it already exists
                if cluster not in data:
                    data[ cluster ] = { "x":[], "y":[] }
                
                # Add entry
                data[ cluster ]["x"].append( temp[0] )
                data[ cluster ]["y"].append( temp[1] )

    # Save it in appropriate files
    for cluster, dataset in data.iteritems():
        with open( outFiles+str(cluster) , "w") as fo:
            for i in range( 0, len(dataset["x"]) ):
                fo.write( dataset["x"][i]+"\t"+dataset["y"][i]+"\n" )
      
    return len(data)
            

# Function for running the actual analysis
def runAnalysis( caseDir ):

    print "Now Doing hier Cluster"
    
    # Get the cluster vs. frame information
    clusterInfo = []
    with open( caseDir+"/analysis/data/cluster_hier_out" , "r") as fi:
        next(fi)
        for line in fi:
            clusterInfo.append( line.split()[1] )
            
    # Separate the dataset into separate files.
    numRmsdDataSets = separateDataSet( 
        caseDir+"/analysis/data/backbone.rmsd", 
        clusterInfo ,
        caseDir+"/analysis/data/cluster_hier_rmsd"
    )
    
    # Create lists of labels and files for plotting
    clusterLabels = []
    clusterFiles = []
    for i in range( 1, numRmsdDataSets+1):
        clusterLabels.append( "Cluster "+str(i) )
        clusterFiles.append( caseDir+"/analysis/data/cluster_hier_rmsd"+str(i) )
    
    # Do the plottin
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "RMSd Cluster Plot", 
        clusterLabels, 
        clusterFiles , 
        "RMSd ($\AA$)",
        scatter = True,
        lowerLegend = True,
        legendLoc = 4
    )