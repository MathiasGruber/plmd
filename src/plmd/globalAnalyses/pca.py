#!/usr/bin/python
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import plmd.generalAnalyses.clusterAnalysis as cluster
import plmd.plotData as myPlot
import os, sys, math
from pylab import plt
import numpy as np
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Function for running the actual analysis
def runAnalysis( caseDirs , resultsDir , noReweight = False):
    
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
            
            ## PCA PLOTTING ON REF DIR PCA COMPONENTS
            #########################################
            
            # Create & run cpptraj for plotting all cases on the axes of the first eigenvector
            # Good URLs for PCA in CPPTRAJ:
            # http://archive.ambermd.org/201404/0243.html
                        
            # PCA plotter
            pcaHandler = pcaFuncs.PCA( 
                resultsDir+"/plots/pcaComparison/PCA_"+caseID+"_on_"+refID+".pdf"
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
            
            ## REWEIGHTING OF PCA PLOTS ON RED DIR PCA COMPONENTS
            #####################################################

            # Check if we should do a reweighted version
            if noReweight == False:
                if os.path.isfile( caseDir+"/md-logs/weights.dat" ):
                    
                    # User info
                    print "aMD weights found. Now attempting 2D reweighting"   
                    
                    # Prepare input file
                    numLines = 0
                    with open(caseDir+"/analysis/data/global_pca", "r") as fi:
                        with open(caseDir+"/analysis/data/global_pca_singleColumn", "w") as fo:
                            next(fi)
                            for line in fi:
                                numLines += 1
                                fo.write( line.split()[1]+"\t"+line.split()[2]+"\n" )

                    # Set the discretization
                    reqBins = 100         
                    discretization = (2*limit) / reqBins                       
                    
                    # Run the reweighting procedure
                    command = "python $PLMDHOME/src/PyReweighting/PyReweighting-2D.py \
                                -input "+caseDir+"/analysis/data/global_pca_singleColumn \
                                -name "+caseDir+"/analysis/data/global_pca_singleColumn_reweighted \
                                -Xdim -"+str(limit)+" "+str(limit)+" \
                                -Ydim -"+str(limit)+" "+str(limit)+" \
                                -discX "+str(discretization)+" \
                                -discY "+str(discretization)+" \
                                -cutoff 1 \
                                -Emax 10 \
                                -job amdweight_CE \
                                -weight "+refDir+"/md-logs/weights.dat | tee -a reweight_variable.log"
                    print "Running command:", command
                    os.system( command )
                    
                    # Create long file for PCA module
                    with open(caseDir+"/analysis/data/global_pca_reweightedDone", "w") as fo:
                        with open(caseDir+"/analysis/data/global_pca_singleColumn_reweighted-pmf-c2.dat", "r") as fi:
                            frame = 0
                            for line in fi:
                                temp = line.split()
                                entries = int(float(temp[2])*10)
                                for i in range(0,entries):
                                    fo.write( str(frame) + "\t" + temp[0] + "\t" + temp[1] +"\n" )
                                    frame += 1

                    # Print block analysis
                    fig, ax = plt.subplots(figsize=(8, 8), nrows=1, ncols=1 )
                    font = {'family' : 'Arial','weight' : 'normal','size' : 10}
                    plt.rc('font', **font)
                    
                    # Now plot the 2d histogram
                    hist = np.load(caseDir+"/analysis/data/global_pca_singleColumn_reweighted_c2EnergyHist.npy")   
                    xedges = np.load(caseDir+"/analysis/data/global_pca_singleColumn_reweighted_c2edgesX.npy")   
                    yedges = np.load(caseDir+"/analysis/data/global_pca_singleColumn_reweighted_c2edgesY.npy")   
                    
                    # Remove points above limit
                    for jy in range(len(hist[0,:])):
                        for jx in range(len(hist[:,0])):
                            if hist[jx,jy] >= 10:
                                hist[jx,jy] = float("inf")
                    
                    # Do plot
                    img = plt.imshow(hist.transpose(),  interpolation='nearest', origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1]] , rasterized=True )
                    
                    # create an axes on the right side of ax. The width of cax will be 5%
                    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)   
                    
                    # Create colorbar
                    colorbar = plt.colorbar(img, ax=ax, cax = cax)
                    colorbar.set_label("Kcal / mol")
                    
                    # Set title, labels etc
                    plt.legend()
                    ax.set_xlabel("PC1", fontsize=12)
                    ax.set_ylabel("PC2", fontsize=12)
                    
                    ax.set_title( "PCA. Case: "+caseID+" Reweighted. Ref case: "+refID )
                    plt.rc('font', **font) 
                    
                    # Save figure
                    fig.savefig(resultsDir+"/plots/pcaComparison/PCA_"+caseID+"_on_"+refID+"_reweighted.pdf")
                    

            ## CLUSTER PLOTS ON PCA COMPONENTS
            ##################################

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
                        clusterType+"_"+caseID+"_on_"+refID, 
                        clusterLabels, 
                        clusterFiles , 
                        "PC2",
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