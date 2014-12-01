#!/usr/bin/python
import plmd.plotData as myPlot
from pylab import plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages
import plmd.generalAnalyses.clusterAnalysis as cluster

# Function for running the actual analysis
def runAnalysis( caseDir , eps, minPoints ):

    # User information
    print "Now Doing dbscan Cluster"
    titlePostpend = "eps: "+str(eps)+", minPoints: "+str(minPoints)
    
    ## KDIST PLOT
    #############

    if os.path.isfile( caseDir+"/analysis/data/Kdist.1.dat" ):
        myPlot.plotData( 
            caseDir+"/analysis/plots" , 
            "Kdist Plots", 
            ["1-dist","2-dist","3-dist","4-dist","5-dist","6-dist"], 
            [caseDir+"/analysis/data/Kdist.1.dat", 
             caseDir+"/analysis/data/Kdist.2.dat", 
             caseDir+"/analysis/data/Kdist.3.dat", 
             caseDir+"/analysis/data/Kdist.4.dat", 
             caseDir+"/analysis/data/Kdist.5.dat", 
             caseDir+"/analysis/data/Kdist.6.dat"] , 
            "k-dist",
            xUnit = "points",
            skipLines = 1,
            xLimits=[0,100])    
    
    # Instantiate the class
    if os.path.isfile(caseDir+"/analysis/data/cluster_dbscan_out"):
        
        # Start the handler
        handler = cluster.clusterBase( caseDir+"/analysis/data/cluster_dbscan_out" )
        
        ## RMSd PLOT
        ############
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
        
        # First one is noise
        clusterLabels[0] = "Noise"   
        
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
        plt.title( "Cluster DBscan Occupancy Fraction, "+titlePostpend )
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
            
            # First one is noise
            clusterLabels[0] = "Noise"        
            
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