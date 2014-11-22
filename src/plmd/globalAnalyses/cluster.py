#!/usr/bin/python
import os, sys, glob, re
from pylab import plt
import numpy as np
import plmd.plotData as myPlot

from MDAnalysis import *
from MDAnalysis.analysis.align import *

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

# Function for running the actual analysis
def runAnalysis( caseDirs, resultsDir, noReweight = False ):
    
    # User info
    print "Doing Cluster analyses."
    showClusters = 10
    clusterTypes = ["hier", "dbscan"]   
    
    # Take each directory as ref directory
    for refDir in caseDirs:
        for refType in clusterTypes:
        
            # Go through each case, except if it's the same
            for caseDir in caseDirs:
                if caseDir != refDir:
                    for caseType in clusterTypes:
                    
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
                                        
                                        # Get RMSd value
                                        refPdb  = Universe( refFile )
                                        casePdb = Universe( caseFile )
                                        alignto(casePdb, refPdb, select="protein and name CA", mass_weighted=True)
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
                                resultsDir+"/clusterComparison" , 
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
                                resultsDir+"/clusterComparison" , 
                                "similars_"+plotName, 
                                resultsDir+"/data/similars_"+plotName+".dat", 
                                "Case "+caseDir.split("/")[-1]+", Cluster ID", 
                                "Case "+refDir.split("/")[-1]+", Cluster ID" , 
                                xColumn=range(1,showClusters+1),
                                yColumn=range(1,showClusters+1)
                            )
        