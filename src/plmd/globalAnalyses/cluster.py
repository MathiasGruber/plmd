#!/usr/bin/python
import os, sys, glob, re
from pylab import plt
import numpy as np
import plmd.plotData as myPlot
from matplotlib.backends.backend_pdf import PdfPages

from MDAnalysis import *
from MDAnalysis.analysis.align import *

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

# Function for running the actual analysis
def runAnalysis( caseDirs, resultsDir, noReweight = False ):
    
    showClusters = 10
    clusterTypes = ["hier", "dbscan"]   

    ## Do GLOBAL CLUSTER ANALYSIS using DBSCAN
    ##########################################

    # User info
    print "Doing global clustering"

    # Run the clustering analysis. It's saved in each case directory, but just use the first one.
    os.system( "$AMBERHOME/bin/cpptraj -p "+caseDirs[0]+"/md-files/peptide_nowat.prmtop -i "+caseDirs[0]+"/ccptraj_analysis_clustering.ptraj" )

    # Move k-dist from dbScan cluster to data directory
    os.system("mv Kdist.*.dat "+resultsDir+"/data/")                                  

    # Plot kdist... man
    myPlot.plotData( 
            resultsDir+"/plots" , 
            "Kdist Plots", 
            ["1-dist","2-dist","3-dist","4-dist","5-dist","6-dist"], 
            [resultsDir+"/data/Kdist.1.dat", 
             resultsDir+"/data/Kdist.2.dat", 
             resultsDir+"/data/Kdist.3.dat", 
             resultsDir+"/data/Kdist.4.dat", 
             resultsDir+"/data/Kdist.5.dat", 
             resultsDir+"/data/Kdist.6.dat"] , 
            "k-dist",
            xUnit = "points",
            skipLines = 1,
            xLimits=[0,100]) 

    ## OCCUPANCY PLOT
    #################
    
    for cluster in clusterTypes:

        # Get the data to plot
        names, fractions = [],[]
        if os.path.isfile( resultsDir+"/data/cluster_"+cluster+"_summary.dat" ):
            with open(resultsDir+"/data/cluster_"+cluster+"_summary.dat","r") as fi:
                next(fi)
                for aline in fi:
                    if aline:
                        values = aline.split()
                        names.append( "Cluster "+values[0] )                      
                        fractions.append( float(values[2]) )
                        
            # Create an array for the interactions
            y_pos = np.arange(len(names))
            
            # Do a bar plot of fractions 'family' : 'Arial',
            pp = PdfPages( resultsDir+"/plots/Cluster_"+cluster+"_occupancy.pdf" )
            font = {
                    'weight' : 'normal',
                    'size'   : 10}    
            fig = plt.figure(figsize=(16,5))
            plt.barh( y_pos, fractions, align = 'center', color = plt.rcParams['axes.color_cycle'][0]  )   
            plt.yticks(y_pos, names)
            ax = fig.gca()
            ax.set_xlabel("Occupied Fraction", fontsize=12)
            ax.set_ylabel("", fontsize=12)
            plt.title( "Cluster "+cluster+" Occupancy Fraction, " )
            plt.rc('font', **font)        
            plt.savefig(pp, format="pdf",dpi=100)
            pp.close()  

    ## Do cluster comparisons between cases
    ##########################################

    # User info
    print "Doing case cluster comparisons & plotting"
    
    # Save information in these arrays
    clusterFiles = []
    
    # Take each directory as ref directory
    for refDir in caseDirs:
        for refType in clusterTypes:
        
            # Go through each case, except if it's the same
            for caseDir in caseDirs:
                if caseDir != refDir:
                    for caseType in clusterTypes:
                    
                        ## CLUSTER COMPARISON 2D PLOTS
                        ##############################                    
                    
                        # The Data Array
                        data = []   
                    
                        # Go through all cluster files for the refDir
                        refCases = 0
                        refFiles = glob.glob( refDir+"/analysis/structures/cluster/cluster_"+refType+"_centroid.c*.pdb" )     
                        for refFile in sorted(refFiles, key=numericalSort):
                            if refCases < showClusters:
                                refCases += 1
                                data.append([])
                                
                                # Get the rmsd to all files in the caseDir
                                caseCases = 0
                                caseFiles = glob.glob( caseDir+"/analysis/structures/cluster/cluster_"+caseType+"_centroid.c*.pdb" )
                                for caseFile in sorted(caseFiles, key=numericalSort):
                                    if caseCases < showClusters:
                                        caseCases += 1
                                        
                                        # Save files
                                        if caseFile not in clusterFiles:
                                            clusterFiles.append(caseFile)
                                        if refFile not in clusterFiles:
                                            clusterFiles.append(refFile)
                                        
                                        # Get RMSd value
                                        refPdb  = Universe( refFile )
                                        casePdb = Universe( caseFile )
                                        alignto(casePdb, refPdb, select="not (name H*)", mass_weighted=True)
                                        rmsdValue    = rmsd(refPdb.atoms.CA.coordinates(), casePdb.atoms.CA.coordinates())
                                        
                                        # Add to data array
                                        data[ refCases-1 ].append( rmsdValue )
                     
                        # Plot if data is not empty
                        if data != []:
                            
                            # Save data to file
                            plotName = refType+":"+refDir.split("/")[-1]+"_vs_"+caseType+":"+caseDir.split("/")[-1]+"_clusterRmsdMap"
                            np.savetxt( resultsDir+"/data/"+plotName+".dat", data, delimiter='\t', fmt='%10.3f')                    
                            
                            # CA rmsd distance map
                            myPlot.plotDataMap( 
                                resultsDir+"/plots/clusterComparison" , 
                                plotName, 
                                resultsDir+"/data/"+plotName+".dat", 
                                "Case "+caseDir.split("/")[-1]+", Cluster ID", 
                                "Case "+refDir.split("/")[-1]+", Cluster ID" , 
                                xColumn=range(1,showClusters+1),
                                yColumn=range(1,showClusters+1)
                            )
                            
                            # Also save a plot showing similar clusters, i.e. RMSd < 0.5
                            data = [[ 1 if j<1 else 0  for j in i] for i in data] 
                            np.savetxt( resultsDir+"/data/similars_"+plotName+".dat", data, delimiter='\t', fmt='%10.3f')                    
                                                        
                            # CA distance map
                            myPlot.plotDataMap( 
                                resultsDir+"/plots/clusterComparison" , 
                                "similars_"+plotName, 
                                resultsDir+"/data/similars_"+plotName+".dat", 
                                "Case "+caseDir.split("/")[-1]+", Cluster ID", 
                                "Case "+refDir.split("/")[-1]+", Cluster ID" , 
                                xColumn=range(1,showClusters+1),
                                yColumn=range(1,showClusters+1)
                            )
                            
    

    
        