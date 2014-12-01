#!/usr/bin/python
import plmd.plotData as myPlot
import numpy as np
import os

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDirs, resultsDir, backbone , noReweight = False ):

    print "Creating plot of dihedrals"

    # Retrieve info on the peptide
    resNames = backbone.resnames()
    
    # Max free energy
    cbMax = 20
    
    # Go through each residue connection
    for i in range( 1, len(resNames) ):
        for angle in ["phi","psi"]:
            
            # User info
            print "Analysing residue "+angle+"_"+str(i)

            # Labels and files
            dataLabels, dataFiles = [],[]
        
            # Go through all directories and set the required files
            for refDir in caseDirs:
                
                # Check for file
                if os.path.isfile( refDir+"/analysis/data/"+angle+"_"+str(i) ):
                    
                    # Load the values
                    values = np.loadtxt( refDir+"/analysis/data/"+angle+"_"+str(i) , usecols=[1])
                    
                    # Create histogram & save in file
                    hist, binEdges = np.histogram( values , bins=60 )
                    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
                    
                    #so that distrib
                    hist=np.add(hist,0.000000000000000001)  
                    
                    #Convert to free energy in Kcal/mol
                    hist=(0.001987*300)*np.log(hist) 
                    
                    # zero value to lowest energy state
                    hist=np.max(hist)-hist  
                    temphist=hist
                    
                    # set infinity free energy values to is cb_max
                    for x in range(len(temphist[:])):
                        if temphist[x] > cbMax:
                            temphist[x]=cbMax   
                            
                    with open(refDir+"/analysis/data/"+angle+"_"+str(i)+"_hist", "w") as fo:
                        for n in range(0, len(temphist)):
                            fo.write( str(bincenters[n])+"\t"+str(temphist[n])+"\n" )
                        
                    # Do plotting
                    dataFiles.append( refDir+"/analysis/data/"+angle+"_"+str(i)+"_hist" )
                    dataLabels.append( refDir.split("/")[-1]+" "+angle+" "+str(i) )
                    
                    # Check if we should do a reweighted version
                    if noReweight == False:
                        if os.path.isfile( refDir+"/md-logs/weights.dat" ):
                            
                            # User info
                            print "aMD weights found. Now attempting reweighting"
                            
                            # Prepare input file
                            numLines = 0
                            with open(refDir+"/analysis/data/"+angle+"_"+str(i), "r") as fi:
                                with open(refDir+"/analysis/data/"+angle+"_"+str(i)+"_singleColumn.rmsd", "w") as fo:
                                    next(fi)
                                    for line in fi:
                                        numLines += 1
                                        fo.write( line.split()[1]+"\n" )
                            
                            # Run the reweighting procedure
                            os.system("python $PLMDHOME/src/PyReweighting/PyReweighting-1D.py \
                                        -input "+refDir+"/analysis/data/"+angle+"_"+str(i)+"_singleColumn.rmsd \
                                        -name "+refDir+"/analysis/data/"+angle+"_"+str(i)+"_reweighted \
                                        -Xdim -180 180 \
                                        -disc 6 \
                                        -cutoff 10 \
                                        -Emax 20 \
                                        -job amdweight_CE \
                                        -weight "+refDir+"/md-logs/weights.dat | tee -a reweight_variable.log")
                                        
                            # Save references for plotting            
                            dataFiles.append( refDir+"/analysis/data/"+angle+"_"+str(i)+"_reweighted-pmf-c2.dat" )
                            dataLabels.append( refDir.split("/")[-1]+" Reweigthed" )
        
            # Do that plotting
            myPlot.plotData( 
                resultsDir+"/plots" , 
                "Dihedrals, "+angle+", ID"+str(i), 
                dataLabels, 
                dataFiles , 
                "E (kcal/mol)" ,
                skipLines = 1
            )
        