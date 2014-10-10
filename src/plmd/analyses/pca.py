#!/usr/bin/python
import os
import numpy as np
import MDAnalysis
from matplotlib.colors import LogNorm
from pylab import plt,hist2d
from matplotlib.backends.backend_pdf import PdfPages

# Get numpy vector, and cut it into a range and return list of min,max bins
def getBinsFromArray( array ):
    arr_max = np.max( array )
    arr_min = np.min( array )
    myTicks = np.arange( arr_min, arr_max, 2 )
    bins = [ [ myTicks[n], myTicks[n+1] ] for n in range(0, len(myTicks)-1 ) ]    
    return bins
    
# Function to check if a value is in a given list of bins
def getBin( bins, value ):
    n = 0
    for bin in bins:
        if float(value) >= bin[0] and float(value) <= bin[1]:
            return n
        n += 1
    return False

# Function for running the actual analysis
def runAnalysis( caseDir, mdTrajectory ):

    # User info
    print "Doing principal component analysis"

    # Output file name and location
    dataDir = caseDir+"/analysis/data/"
    plotDir = caseDir+"/analysis/plots/"
    structDir = caseDir+"/analysis/structures/"

    # Do the four principal vectors
    frames = []
    pc1 = []
    pc2 = []
    pc3 = []
    pc4 = []

    # Go through analysis file and get eigenvalues
    eigenValues = []
    eigenValueTotal = 0
    with open(dataDir+"evecs-ca.dat","r") as f:
        lookForVec = 1
        for line in f:
            temp = line.split()
            if len(temp) == 2 and temp[0] == str(lookForVec):
                lookForVec += 1
                eigenValueTotal += float(temp[1])
                eigenValues.append(float(temp[1]))
    eigenValues = (np.array(eigenValues) / eigenValueTotal) * 100               
                

    # Get the file with all the projection data
    pcaFile = open(dataDir+"pca12-ca","r")
    n = 0
    for aline in pcaFile:
        if n > 1 and aline:
            values = aline.split()
            frames.append( int(values[0]) )
            pc1.append( float(values[1]) )
            pc2.append( float(values[2]) )
            pc3.append( float(values[3]) )
            pc4.append( float(values[4]) )
        n = n + 1

    # Create numpy arrays
    frames = np.array( frames )
    np_arrays = []
    np_arrays.append( np.array( pc1 ) )
    np_arrays.append( np.array( pc2 ) )
    np_arrays.append( np.array( pc3 ) )
    np_arrays.append( np.array( pc4 ) )
    
    # Do the plotting
    pp = PdfPages( plotDir+"PCA_analysis.pdf" )
        
    # Set the plotting font and default size
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 10}
    
    # Do a plot for each PCA
    for component in range( 1, len(np_arrays) ):

        # User Info        
        print "Plotting component 1 vs. "+str(component)
    
        # Normalize the data DeltaG = -kb T * [ ln( P(v1,v2) ) - ln Pmax ]
        boltzman = 0.0019872041
        temperature = 300
        
        # Plot both distribution & Energy Landscape
        for plotType in [ "energy", "distribution" ]:        
        
            # Create figure
            fig, ax = plt.subplots()                

            # Do the plotting
            if plotType == "energy":
                
                # Create the histogram without plotting, so we can set the units properly        
                H, xedges, yedges = np.histogram2d(np_arrays[component], np_arrays[0], bins=100 )
                H_normalized = H/len(np_arrays[0])
                H = -1 * boltzman * temperature * (np.log( H_normalized )-np.log(np.max(H_normalized)))
                
                # Now plot the 2d histogram
                plt.imshow(H,  interpolation='nearest', origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1]])
                colorbar = plt.colorbar()
                colorbar.set_label("Kcal / mol")
                
            elif plotType == "distribution":
        
                # Directly do the 2d histogram of matplotlib        
                hist2d(np_arrays[0], np_arrays[component], bins=100)
                colorbar = plt.colorbar()
                colorbar.set_label("Occurances")
        
            # Set title, labels etc
            plt.legend()
            ax = fig.gca()
            ax.set_xlabel("PC1 ({0:.2f}%)".format(eigenValues[0]), fontsize=12)
            ax.set_ylabel("PC"+str(component+1)+" ({0:.2f}%)".format(eigenValues[component]), fontsize=12)
            plt.title( "PCA Analysis" )
            plt.rc('font', **font)   
    
            # Save figure in pdf and png
            plt.savefig(pp, format="pdf",dpi=300)
        
   
    
    # Close the figure
    pp.close()
    
    # Go through the vectors X, and for vector1, vectorX,
    # find a frame representing the given point and output pdb files
    ################################################################
    
    # Get the bins for the first component, which is always plotted
    primaryBins = getBinsFromArray( np_arrays[0] )           
    
    # Select everything for printing in pdb file
    selection = mdTrajectory.selectAtoms("all")    
    
    # Go through other components
    for component in range( 1,len( np_arrays ) ):
        
        # Get the bins for the current component
        secondaryBins = getBinsFromArray( np_arrays[ component ] )
        
        # 2d bin matrix
        binMatrix = np.ones( ( len(primaryBins), len(secondaryBins) ) )
        
        # Loop through file and find 
        pcaFile = open(dataDir+"pca12-ca","r")
        n = 0
        for aline in pcaFile:
            if n > 1 and aline:
                values = aline.split()
                
                # Check the frame is in the first component
                primaryBin = getBin( primaryBins, values[1] )
                if primaryBin != False:
                    
                    # Check that the frame is in the secondary component bin
                    secondaryBin = getBin( secondaryBins, values[component+1] )
                    if secondaryBin != False:
                        
                        # Check if this bin was already analysed
                        if binMatrix[ primaryBin, secondaryBin ] == 1:

                            # Write to ptraj buffer
                            bin1 = (primaryBins[ primaryBin ][0]+primaryBins[ primaryBin ][1])*0.5
                            bin2 = (secondaryBins[ secondaryBin ][0]+secondaryBins[ secondaryBin ][1])*0.5
                            
                            # Use MDAnalysis to pring PDB
                            filename = structDir+"pca"+str(component+1)+"/PC1_"+str(bin1)+"_PC"+str(component+1)+"_"+str(bin2)+"_Frame_"+values[0]+".pdb"
                            
                            # Print stuff
                            print "Looking at frame: ", values[0]
                            print "Writing to file: ", filename
                            print "Current timestep: ", mdTrajectory.trajectory[ int(values[0]) ]
                                
                            # Write the PDB file
                            W = MDAnalysis.Writer( filename, remarks='Have a nice day', )
                            W.write(selection)
                            W.close()
                            
                            # Open it and remove REMARK
                            lines = []
                            with open( filename , "r" ) as fl:
                                lines = fl.readlines()
                            with open( filename , "w" ) as fl:
                                for line in lines:
                                    if "REMARK" not in line:
                                        fl.write(line)
                            
                            # These are assigned, so remove the bins
                            binMatrix[ primaryBin, secondaryBin ] = 0
                        
            n = n + 1
    


