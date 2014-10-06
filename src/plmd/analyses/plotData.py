#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator

# Set plotting colors
# plt.rcParams['axes.color_cycle'] = ['darkred', 'orange','green', 'yellow', 'lightgreen',  'lightblue']

# Plot graphs
# Can plot multiple files
def plotData( caseDir , title, labels , files , unit , xUnit = "Time (ps)", types=None, skipLines = 0, xFactor = 1 ):
    
    # Output file name and location
    dataDir = caseDir+"/analysis/data/"
    plotDir = caseDir+"/analysis/plots/"
    
    # plot using pdf files
    pp = PdfPages( plotDir+title+".pdf" )
    
    # Set figure size
    fig = plt.figure(figsize=(8,6))
    
    # Go through the files to plot
    i = 0
    for filename in files:
        xData, yData = [],[]
        qbfile = open(dataDir+filename,"r")
        n = 0
        for aline in qbfile:
            if n > skipLines:
                values = aline.split()
                xData.append(float(values[0])*xFactor)
                yData.append(values[1])
            n = n + 1

        # Set linetype for the plot
        lineType = '-'
        if types != None:
            lineType = types[i]
            
        # Do the plotting (rasterize to reduce load time)
        plt.plot( xData, yData , lineType, label = labels[i] , rasterized=True)
        
        # Next in line
        i = i + 1
                
    # Set title, labels etc
    plt.legend()
    ax = fig.gca()
    ax.set_xlabel(xUnit, fontsize=12)
    ax.set_ylabel(unit, fontsize=12)
    plt.title( title )

    # Set the plotting font and default size
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 10}
    plt.rc('font', **font)        
    
    # Save figure in pdf and png
    plt.savefig(pp, format="pdf",dpi=300)
    
    # Close the figure
    pp.close()

# Do 2d plots
def plotDataMap( caseDir , title , file , xUnit , yUnit , xColumn=0 , yColumn=0 , skipRow = 0, skipColumn = 0 ,xFactor=1, yFactor=1):

    # Output file name and location
    dataDir = caseDir+"/analysis/data/"
    plotDir = caseDir+"/analysis/plots/"
    
    # plot using pdf files
    pp = PdfPages( plotDir+title+".pdf" )
    
    # Create figure
    fig, ax = plt.subplots()
    
    # Load data
    data = []
    qbfile = open(dataDir+file,"r")
    n = 0
    for aline in qbfile:
        if n != (skipRow-1):
            strings = aline.split()
            if skipColumn != 0:
                strings.pop((skipColumn-1))
            numbers = map(float, strings)
            data.append( numbers )
        n = n + 1
    data = np.array(data)

    # Create heatmap in blue
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues , rasterized=True)

    # Factor x axis
    labels = ax.get_xticks().tolist()
    i = 0
    for label in labels:
        labels[i] = str( float(labels[i])*xFactor )
        i = i + 1
    ax.set_xticklabels(labels)

    # Factor y axis
    labels = ax.get_yticks().tolist()
    i = 0
    for label in labels:
        labels[i] = str( float(labels[i])*yFactor )
        i = i + 1
    ax.set_yticklabels(labels)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    # Set it up
    if xColumn != 0 and yColumn != 0:
        ax.set_xticklabels(xColumn, minor=False)
        ax.set_yticklabels(yColumn, minor=False)
    else:
        ax.xaxis.set_major_locator(LinearLocator(numticks=len(labels)))
        ax.yaxis.set_major_locator(LinearLocator(numticks=len(labels)))

    # Set title
    ax.set_xlabel(xUnit, fontsize=12)
    ax.set_ylabel(yUnit, fontsize=12)
    myBar = plt.colorbar(heatmap)
    myBar.set_label("Distance ($\AA$)")
    plt.title( title )
    
    # Set the plotting font and default size
    font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
    plt.rc('font', **font)

    # Save figure
    plt.savefig(pp, format="pdf",dpi=300)

    # Close the figure
    pp.close()
