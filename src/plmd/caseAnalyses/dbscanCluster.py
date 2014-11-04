#!/usr/bin/python
import plmd.plotData as myPlot
from pylab import plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import plmd.generalAnalyses.clusterAnalysis as cluster

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # User information
    print "Now Doing dbscan Cluster"
    
    # Instantiate the class
    handler = cluster.clusterBase( caseDir+"/analysis/data/cluster_dbscan_out" )

    # RMSd PLOT
    numRmsdDataSets = handler.separateDataSet( 
        caseDir+"/analysis/data/backbone.rmsd",
        caseDir+"/analysis/data/cluster_dbscan_rmsd_"
    )
    
    # User info
    print "Number of datasets for RMSd: "+str(numRmsdDataSets) 
    
    # Create lists of labels and files for plotting
    clusterLabels = []
    clusterFiles = []
    for i in range( 0, numRmsdDataSets):
        clusterLabels.append( "Cluster "+str(i) )
        clusterFiles.append( caseDir+"/analysis/data/cluster_dbscan_rmsd_d1_c"+str(i) )
    
    # Do the plottin
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "RMSd DBscan Cluster Plot", 
        clusterLabels, 
        clusterFiles , 
        "RMSd ($\AA$)",
        scatter = True,
        legendLoc = 4
    )
    
    ## OCCUPANCY PLOT
    #################
    
    # Get the data to plot
    names, fractions = [],[]
    with open(caseDir+"/analysis/data/cluster_dbscan_summary.dat","r") as fi:
        next(fi)
        for aline in fi:
            if aline:
                values = aline.split()
                names.append( "Cluster "+values[0] )                      
                fractions.append( float(values[2]) )
                
    # Create an array for the interactions
    y_pos = np.arange(len(names))
    
    # Do a bar plot of fractions
    pp = PdfPages( caseDir+"/analysis/plots/ClusterDBscanOccupancy.pdf" )
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 10}    
    fig = plt.figure(figsize=(16,5))
    plt.barh( y_pos, fractions, align = 'center', color = plt.rcParams['axes.color_cycle'][0]  )   
    plt.yticks(y_pos, names)
    ax = fig.gca()
    ax.set_xlabel("Occupied Fraction", fontsize=12)
    ax.set_ylabel("", fontsize=12)
    plt.title( "Cluster DBscan Occupancy Fraction" )
    plt.rc('font', **font)        
    plt.savefig(pp, format="pdf",dpi=100)
    pp.close()    
    
    ## PCA PLOT
    ###########
    
    # Separate the dataset
    numPCAdataSets = handler.separateDataSet( 
        caseDir+"/analysis/data/pca12-ca", 
        caseDir+"/analysis/data/cluster_dbscan_pca_",
        xColumn = 1
    )    
    
    # User info
    print "Number of datasets for PCA: "+str(numPCAdataSets)    
    
    # Go through each PCA component
    for n in range(2,4):    
    
        # Create lists of labels and files for plotting
        clusterLabels = []
        clusterFiles = []
        for i in range( 0, numPCAdataSets):
            clusterLabels.append( "Cluster "+str(i) )
            clusterFiles.append( caseDir+"/analysis/data/cluster_dbscan_pca_d"+str(n)+"_c"+str(i) )
        
        print "Calling plot"
        print "labels", clusterLabels  
        print "files", clusterFiles  
        
        myPlot.plotData( 
            caseDir+"/analysis/plots" , 
            "PCA DBscan Cluster, 1vs"+str(n), 
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