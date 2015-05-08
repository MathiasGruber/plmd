#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LinearLocator

# Plot graphs
# Can plot multiple files
def plotData( 
    outputDir , 
    title, 
    labels , 
    inputFiles , 
    unit , 
    xUnit = "Time (ps)", 
    types = None, 
    markTypes = None, 
    scatter = False, 
    tightXlimits = True,
    tightYlimits = False,
    skipLines = 0, 
    xFactor = 1,
    showTitle = True,
    figWidth = 8,
    figHeight = 6,
    legendLoc = 1,
    legendAlpha = 0,
    legendFrame = 0,
    xLimits = False,
    yLimits = False,
    logY = False,
    fontSize = 10
):
    
    # plot using pdf 
    pp = PdfPages( outputDir+"/"+title+".pdf" )
    
    # Set figure size
    fig = plt.figure(
        figsize=(figWidth,figHeight)
    )
    
    # Go through the to plot
    i = 0
    xmin,xmax,ymin,ymax = 0,0,0,0
    
    for filename in inputFiles:
        xData, yData = [],[]
        qbfile = open(filename,"r")
        n = 0
        for aline in qbfile:
            if n >= skipLines:
                values = aline.split()
                if yLimits == False or float(values[1]) > yLimits[0]:
                    xData.append(float(values[0])*xFactor)
                    yData.append(float(values[1]))
            n = n + 1
            
        # Check if there's any data
        if yData:            
        
            # Limit data to 500 data points
            if n > 1000:
                timesOver = int(math.floor( n / 1000. ))
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
            
                # Set linetype for the plot
                markType = 'o'
                if markTypes != None:
                    markType = markTypes[i]
            
                # Do the plot
                colors = plt.rcParams['axes.color_cycle']
                color = colors[i % len(colors)]
                plt.scatter( xData, yData , s=8, color=color, label = labels[i] , rasterized=True, marker = markType)
            
            # Get limits
            xmax = np.max( xData ) if np.max( xData ) > xmax else xmax
            xmin = np.min( xData ) if np.min( xData ) < xmin else xmin
            ymax = np.max( yData ) if np.max( yData ) > ymax else ymax
            ymin = np.min( yData ) if np.min( yData ) < ymin else ymin
    
        # Next in line
        i = i + 1
        
    # Set title, labels etc
    ax = fig.gca()
    ax.set_xlabel(xUnit, fontsize=fontSize+2)
    ax.set_ylabel(unit, fontsize=fontSize+2)
    
    # Check if tight x-limits
    if tightXlimits == True:
        xmin,xmax = math.floor(xmin), math.ceil(xmax)
        print "Now setting X limites to: ",xmin, xmax
        ax.set_xlim([xmin,xmax]) 
        
    # Check if tight y-limits
    if tightYlimits == True:
        ymin,ymax = math.floor(ymin), math.ceil(ymax)
        print "Now setting Y limites to: ",ymin, ymax
        ax.set_ylim([ymin,ymax]) 
        
    # Force limits
    if xLimits != False:
        ax.set_xlim(xLimits) 
    if yLimits != False:
        ax.set_ylim(yLimits) 
        
    # Log scale
    if logY != False:
        ax.set_yscale('log')
        
    # Plot title
    if showTitle == True:
        plt.title( title , fontsize=fontSize+4 )
    
    # Create the legend
    plt.legend(
        loc=legendLoc, 
        ncol=1, 
        framealpha=legendAlpha, 
        frameon=legendFrame
    )   
    
    # Set the plotting font and default size 'family' : 'Arial',
    font = {
            'weight' : 'normal',
            'size'   : fontSize}
    plt.rc('font', **font)        
    
    # Save figure in pdf and png
    pp.savefig(
        fig,
        dpi=150
    )
    
    # Close the figure
    pp.close()

# Do 2d plots
def plotDataMap( outputDir , title , inputFile , xUnit , yUnit , xColumn=0 , yColumn=0 , skipRow = 0, skipColumn = 0 ,xFactor=1, yFactor=1, fontSize=10 ):

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
    ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

    # Set it up
    if xColumn != 0 and yColumn != 0:
        ax.set_xticklabels(xColumn, minor=False)
        ax.set_yticklabels(yColumn, minor=False)
    else:
        ax.xaxis.set_major_locator(LinearLocator(numticks=len(labels)))
        ax.yaxis.set_major_locator(LinearLocator(numticks=len(labels)))

    # Set title
    ax.set_xlabel(xUnit, fontsize=fontSize+2)
    ax.set_ylabel(yUnit, fontsize=fontSize+2)
    myBar = plt.colorbar(heatmap)
    myBar.set_label("Distance ($\AA$)")
    plt.title( title  , fontsize=fontSize+4 )
    
    # Set the plotting font and default size 'family' : 'Arial',
    font = {
        'weight' : 'normal',
        'size'   : fontSize}
    plt.rc('font', **font)

    # Save figure
    plt.savefig(pp, format="pdf",dpi=150)

    # Close the figure
    pp.close()
