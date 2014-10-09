#!/usr/bin/python
import numpy as np
from matplotlib.colors import LogNorm
from pylab import plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # User info
    print "Doing principal component analysis"

    # Output file name and location
    dataDir = caseDir+"/analysis/data/"
    plotDir = caseDir+"/analysis/plots/"

    # Do the four principal vectors
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
                

    # Get the file with all the data
    pcaFile = open(dataDir+"pca12-ca","r")
    n = 0
    for aline in pcaFile:
        if n > 1 and aline:
            values = aline.split()
            pc1.append( float(values[1]) )
            pc2.append( float(values[2]) )
            pc3.append( float(values[3]) )
            pc4.append( float(values[4]) )
        n = n + 1

    # Create numpy arrays
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
        
        # Create figure
        fig, ax = plt.subplots()        

        # Create the histogram without plotting, so we can set the units properly        
        H, xedges, yedges = np.histogram2d(np_arrays[0], np_arrays[component], bins=100 )
        H_normalized = H/len(np_arrays[0])
        H = -1 * boltzman * temperature * (np.log( H_normalized )-np.log(np.max(H_normalized)))
        
        # Now plot the 2d histogram
        plt.imshow(H,  interpolation='nearest')
        colorbar = plt.colorbar()
        colorbar.set_label("Kcal / mol")
        
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
    
    


