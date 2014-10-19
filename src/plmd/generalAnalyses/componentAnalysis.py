#!/usr/bin/python
import math
import numpy as np
import MDAnalysis
from pylab import plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

class PCA():
    
    # Constructor
    def __init__( self , outputFile, subplotRows = 1, subplotColums = 1 ):
        
        # File for plotting
        self.pp = PdfPages( outputFile )
        
        # Save the rows and columns
        self.rows, self.columns = int(subplotRows), int(subplotColums)
        self.spots = self.rows*self.columns   
        
        # The main data array. Used also for saving structures
        self.np_arrays =  []
        
        # Counter for the amount of figures plotted
        self.subPlots = 0
        
        # Array for storing the colobars of the plots (these need to be removed)
        self.colorBars = []
        
    # Get active ax.
    def getActiveAx( self ):
        if self.rows == 1 and self.columns == 1:
            return self.axarr
        elif self.rows == 1 or self.columns == 1:
            return self.axarr[ self.subPlots ]
        else:
            # Determine row and column of this plot
            pickedRow = math.floor( self.subPlots / self.columns )
            pickedColumn = self.subPlots - pickedRow*self.columns
            return self.axarr[ pickedRow ][ pickedColumn ]
    
    # Close the figure
    def savePlot( self ):

        # Save any hanging plots
        if self.subPlots > 0:
            self.savePdfPage()      
        
        # The end
        self.pp.close()
        
    # Create subplot page
    def createPdfPage( self ):
        print "Creating a new subplot grid"
        self.fig, self.axarr = plt.subplots( nrows=self.rows, ncols=self.columns )
        
    # Save current plot
    def savePdfPage( self ):

        # Delete any unused subplots        
        if self.subPlots < self.spots:
            print self.axarr
            for i in range( self.subPlots, self.spots ):
                print "Deleting subplot with id: ",i
                self.fig.delaxes( self.getActiveAx() )
                self.subPlots += 1
        
        # Ensure tight layout so everything is visible
        plt.tight_layout()
        
        # Do the save
        print "Saving pdf page"
        self.pp.savefig( self.fig, dpi=100 )
        
        # Clear figures
        plt.cla()
        
        # We have to explicitly clear color bars due to the way they are added
        self.clearColorbars()
    
    def clearColorbars(self):
        if len(self.colorBars) > 0:
            for bar in self.colorBars:
                bar.remove()
            self.colorBars = []
        

    # Function for running the actual analysis
    def plotPCA( self, plotTitleIdentifier, dataDir, eigenVectorFile , eigenValueFile = False , eigenVectorCount = 4, plotDistibution = True ):
    
        # If this is the first plot, create grid
        if self.subPlots == 0:
            self.createPdfPage()
    
        # User info
        print "Doing principal component analysis"
    
        # Do the four principal vectors
        frames = []

        # Principal components
        pcs = []
        for i in range( 0, eigenVectorCount ):
            pcs.append([])
    
        # Go through analysis file and get eigenvalues
        if eigenValueFile != False:
            eigenValues = []
            eigenValueTotal = 0
            with open(dataDir+eigenValueFile,"r") as f:
                lookForVec = 1
                for line in f:
                    temp = line.split()
                    if len(temp) == 2 and temp[0] == str(lookForVec):
                        lookForVec += 1
                        eigenValueTotal += float(temp[1])
                        eigenValues.append(float(temp[1]))
            eigenValues = (np.array(eigenValues) / eigenValueTotal) * 100               
    
        # Get the file with all the projection data
        pcaFile = open(dataDir+eigenVectorFile,"r")
        n = 0
        for aline in pcaFile:
            if n > 1 and aline:
                values = aline.split()
                
                # Add frames
                frames.append( int(values[0]) )
                
                # If requesting more vectors than present
                if eigenVectorCount > len(values):
                    raise Exception("CPPTRAJ has not projected "+str(eigenVectorCount)+" vectors")
                
                # Add to vectors
                for i in range( 0, eigenVectorCount ):
                    pcs[i].append( float(values[i+1]) )
            n = n + 1
    
        # Create numpy arrays
        frames = np.array( frames )
        self.np_arrays = [ np.array( pc ) for pc in pcs ]
        
        # Set the plotting font and default size
        font = {'family' : 'Arial',
                'weight' : 'normal',
                'size'   : 10}
                
        # Do a plot for each PCA
        for component in range( 1, len(self.np_arrays) ):
            
            # User Info        
            print "Plotting component 1 vs. "+str(component)
        
            # Normalize the data DeltaG = -kb T * [ ln( P(v1,v2) ) - ln Pmax ]
            boltzman = 0.0019872041
            temperature = 300
            
            # Plot both distribution & Energy Landscape
            for plotType in [ "energy", "distribution" ] if plotDistibution else [ "energy" ]:        
            
                # New subplot
                ax = self.getActiveAx()
                
                # Increase subplot counter
                self.subPlots += 1
    
                # Do the plotting
                if plotType == "energy":
                    
                    # Create the histogram without plotting, so we can set the units properly        
                    H, xedges, yedges = np.histogram2d(self.np_arrays[component], self.np_arrays[0], bins=100 )
                    H_normalized = H/len(self.np_arrays[0])
                    H = -1 * boltzman * temperature * (np.log( H_normalized )-np.log(np.max(H_normalized)))
                    
                    # Now plot the 2d histogram
                    img = ax.imshow(H,  interpolation='nearest', origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1]] , rasterized=True )
                    colorbar = plt.colorbar(img, ax=ax)
                    colorbar.set_label("Kcal / mol")
                    self.colorBars.append(colorbar)
                    
                elif plotType == "distribution":
            
                    # Directly do the 2d histogram of matplotlib        
                    _, _, _, img = ax.hist2d(self.np_arrays[0], self.np_arrays[component], bins=100 , rasterized=True, norm=LogNorm() )
                    colorbar = plt.colorbar(img, ax=ax)
                    colorbar.set_label("Occurances")
                    self.colorBars.append(colorbar)
            
                mini = np.abs(np.min( [np.min(self.np_arrays[0]), np.min(self.np_arrays[component])] ))  
                maxi = np.abs(np.max( [np.max(self.np_arrays[0]), np.max(self.np_arrays[component])] ))
                limits = int(math.ceil(np.max( [mini,maxi] )))
                print limits, "BLAH BLAH BLAH"                
                
                plt.ylim([-limits,limits])
                plt.xlim([-limits,limits])            
            
                # Set title, labels etc
                plt.legend()
                if eigenValueFile != False:
                    ax.set_xlabel("PC1 ({0:.2f}%)".format(eigenValues[0]), fontsize=12)
                    ax.set_ylabel("PC"+str(component+1)+" ({0:.2f}%)".format(eigenValues[component]), fontsize=12)
                else:
                    ax.set_xlabel("PC1", fontsize=12)
                    ax.set_ylabel("PC"+str(component+1), fontsize=12)
                
                ax.set_title( "PCA Analysis. "+plotTitleIdentifier )
                #ax.rc('font', **font)   
        
                # Save pdf page if it's filled
                if self.subPlots >= (self.rows*self.columns):
                    print "Now saving to PDF. Number of plots: ",self.subPlots
                    self.savePdfPage()
                    self.subPlots = 0
        
        
    # Go through the vectors X, and for vector1, vectorX,
    # find a frame representing the given point and output pdb files
    ################################################################
    def saveStructures( self, mdTrajectory, structDir, eigenVectorFile ):
        
        # Get the bins for the first component, which is always plotted
        primaryBins = self.getBinsFromArray( self.np_arrays[0] )           
        
        # Select everything for printing in pdb file
        selection = mdTrajectory.selectAtoms("all")    
        
        # Go through other components
        for component in range( 1,len( self.np_arrays ) ):
            
            # Get the bins for the current component
            secondaryBins = self.getBinsFromArray( self.np_arrays[ component ] )
            
            # 2d bin matrix
            binMatrix = np.ones( ( len(primaryBins), len(secondaryBins) ) )
            
            # Loop through file and find 
            pcaFile = open(eigenVectorFile,"r")
            n = 0
            for aline in pcaFile:
                if n > 1 and aline:
                    values = aline.split()
                    
                    # Check the frame is in the first component
                    primaryBin = self.getBin( primaryBins, values[1] )
                    if primaryBin != False:
                        
                        # Check that the frame is in the secondary component bin
                        secondaryBin = self.getBin( secondaryBins, values[component+1] )
                        if secondaryBin != False:
                            
                            # Check if this bin was already analysed
                            if binMatrix[ primaryBin, secondaryBin ] == 1:
    
                                # Write to ptraj buffer
                                bin1 = (primaryBins[ primaryBin ][0]+primaryBins[ primaryBin ][1])*0.5
                                bin2 = (secondaryBins[ secondaryBin ][0]+secondaryBins[ secondaryBin ][1])*0.5
                                
                                # Use MDAnalysis to pring PDB
                                filename = structDir+"pca"+str(component+1)+"/PC1_"+str(bin1)+"_PC"+str(component+1)+"_"+str(bin2)+"_Frame_"+values[0]+".pdb"
                                
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
        
    # Get numpy vector, and cut it into a range and return list of min,max bins
    def getBinsFromArray( self, array ):
        arr_max = np.max( array )
        arr_min = np.min( array )
        myTicks = np.arange( arr_min, arr_max, 2 )
        bins = [ [ myTicks[n], myTicks[n+1] ] for n in range(0, len(myTicks)-1 ) ]    
        return bins
        
    # Function to check if a value is in a given list of bins
    def getBin( self, bins, value ):
        n = 0
        for bin in bins:
            if float(value) >= bin[0] and float(value) <= bin[1]:
                return n
            n += 1
        return False


