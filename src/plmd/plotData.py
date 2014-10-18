#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator


# Set plotting colors
#plt.rcParams['axes.color_cycle'] = ['darkred', 'orange','green', 'yellow', 'lightgreen',  'lightblue']
plt.rcParams['axes.color_cycle'] = ["348ABD", "7A68A6", "A60628", "467821", "CF4457", "188487", "E24A33"]

# Plot graphs
# Can plot multiple files
def plotData( outputDir , title, labels , inputFiles , unit , xUnit = "Time (ps)", types=None, scatter = False, skipLines = 0, xFactor = 1 ):
    
    # plot using pdf 
    pp = PdfPages( outputDir+"/"+title+".pdf" )
    
    # Set figure size
    fig = plt.figure(figsize=(8,6))
    
    # Go through the to plot
    i = 0
    for filename in inputFiles:
        xData, yData = [],[]
        qbfile = open(filename,"r")
        n = 0
        for aline in qbfile:
            if n > skipLines:
                values = aline.split()
                xData.append(float(values[0])*xFactor)
                yData.append(values[1])
            n = n + 1
            
        # Limit data to 500 data points
        if n > 500:
            timesOver = int(math.floor( n / 500. ))
            if timesOver > 2:
                xData = xData[::timesOver]
                yData = yData[::timesOver]
        
        # If xFactor is above 1, save the corrected data set in a new data file
        if xFactor > 1:
            with open( filename+"_timeCorrected", "w" ) as fo:
                for n in range( 0, len(xData) ):
                    fo.write( str(xData[n]) + "\t" + str(yData[n])  +"\n")
        
        # Do the plotting (rasterize to reduce load time)
        if scatter == False:
            
             # Set linetype for the plot
            lineType = '-'
            if types != None:
                lineType = types[i]
            
            # Do the plot
            plt.plot( xData, yData , lineType, label = labels[i] , rasterized=True)
        
        else:
        
            # Do the plot
            color = plt.rcParams['axes.color_cycle'][i]
            area = 2
            plt.scatter( xData, yData , s=area, color=color, label = labels[i] , rasterized=True)
        
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
    plt.savefig(pp, format="pdf",dpi=150)
    
    # Close the figure
    pp.close()

# Do 2d plots
def plotDataMap( outputDir , title , inputFile , xUnit , yUnit , xColumn=0 , yColumn=0 , skipRow = 0, skipColumn = 0 ,xFactor=1, yFactor=1):

    # plot using pdf files
    pp = PdfPages( outputDir+"/"+title+".pdf" )
    
    # Create figure
    fig, ax = plt.subplots()
    
    # Load data
    data = []
    qbfile = open(inputFile,"r")
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
    plt.savefig(pp, format="pdf",dpi=150)

    # Close the figure
    pp.close()
