#!/usr/bin/python
import numpy as np
from pylab import plt
from matplotlib.backends.backend_pdf import PdfPages

# Function for running the actual analysis
def runAnalysis( caseDir ):

    # User info 
    print "Plotting hydrogren bonds"

    # Output file name and location
    dataDir = caseDir+"/analysis/data/"
    plotDir = caseDir+"/analysis/plots/"

    # Do a histogram plot of Hbond information
    pp = PdfPages( plotDir+"hbonds.pdf" )
    
    # Get the data to plot
    names, fractions, avgDists, avgAngles = [],[],[],[]
    qbfile = open(dataDir+"hbond.avg","r")
    n = 0
    for aline in qbfile:
        if n > 1 and n < 17:
            if aline:
                values = aline.split()
                names.append( values[0]+" - "+values[1] )                      
                fractions.append( float(values[4]) )
                avgDists.append( float(values[5]) )
                avgAngles.append( float(values[6]) )
        n = n + 1

    # Create an array for the interactions
    y_pos = np.arange(len(names))
    
    # Set the plotting font and default size
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 10}
        
    # Do a bar plot of fractions
    fig = plt.figure(figsize=(16,5))
    plt.barh(y_pos, fractions, align='center')   
    plt.yticks(y_pos, names)
    ax = fig.gca()
    ax.set_xlabel("Occupied Fraction", fontsize=12)
    ax.set_ylabel("H-bond Interaction", fontsize=12)
    plt.title( "Hydrogen Bonds: Occupancy" )
    plt.rc('font', **font)        
    plt.savefig(pp, format="pdf",dpi=300)

    # Do a bar plot of avg dists
    fig = plt.figure(figsize=(16,5))
    plt.barh(y_pos, avgDists, align='center')   
    plt.yticks(y_pos, names)
    ax = fig.gca()
    ax.set_xlabel("Average Distance", fontsize=12)
    ax.set_ylabel("H-bond Interaction", fontsize=12)
    plt.title( "Hydrogen Bonds: Average Distance" )
    plt.rc('font', **font)        
    plt.savefig(pp, format="pdf",dpi=300)
    
    # Do a bar plot of avg angles
    fig = plt.figure(figsize=(16,5))
    plt.barh(y_pos, avgAngles, align='center')   
    plt.yticks(y_pos, names)
    ax = fig.gca()
    ax.set_xlabel("Average Angle", fontsize=12)
    ax.set_ylabel("H-bond Interaction", fontsize=12)
    plt.title( "Hydrogen Bonds: Average Angle" )
    plt.rc('font', **font)        
    plt.savefig(pp, format="pdf",dpi=300)
    
    # Close the figure
    pp.close()
