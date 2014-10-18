#!/usr/bin/python
import MDAnalysis
import numpy as np
import plmd.plotData as myPlot
from pylab import plt,hist2d
from matplotlib.backends.backend_pdf import PdfPages

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , backbone , timeFactor ):

    # Retrieve info on the peptide
    resNames = backbone.resnames()
    
    # Go through each residue connection
    for i in range( 1, len(resNames) ):
        
        # User info
        print "Plotting dihedrals for residue: "+str(i)
        
        # Paths for the two files
        psiPath = caseDir+"/analysis/data/psi_"+str(i)  
        phiPath = caseDir+"/analysis/data/phi_"+str(i)

        # Common Plot command
        myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "Dihedral Angles. Res ID: "+str(i)+", Residue name: "+str(resNames[i]), 
        ["$\Psi_"+str(i)+"$","$\Phi_"+str(i)+"$"],
        [ psiPath , phiPath ] , 
        "Angle (Degrees)", 
        xFactor = timeFactor )
 
        # Create a Ramachandran plot
        ############################

        # User info
        print "Creating a Ramachandran plot for residue: "+str(i) 
 
        # Get the components
        components = []
        for path in [ phiPath, psiPath ]:
            components.append([])
            index = len(components)-1
            with open(path, "r") as fi:
                lines = fi.readlines()
                for line in lines:
                    temp = line.split()
                    try:
                        components[ index ].append( float( temp[1] ) )
                    except ValueError:
                        print "Reading Header: ",line
        
        # Set to numpy
        np_arrays = [ np.array( component ) for component in components ]  
        
        # Do the plotting
        title = "Ramachandran Plot. Res ID: "+str(i)+", Residue name: "+str(resNames[i])
        pp = PdfPages( caseDir+"/analysis/plots/"+title+".pdf" )
        fig = plt.figure() #figsize=(8,6)
        ax = fig.gca()
        ax.set_xlabel("$\Phi_"+str(i)+"$", fontsize=12)
        ax.set_ylabel("$\Psi_"+str(i)+"$", fontsize=12)
        
        # Create the histogram without plotting, so we can set the units properly   
        boltzman = 0.0019872041
        temperature = 300
        H, xedges, yedges = np.histogram2d(np_arrays[1], np_arrays[0], bins=100 )
        H_normalized = H/len(np_arrays[0])
        H = -1 * boltzman * temperature * (np.log( H_normalized )-np.log(np.max(H_normalized)))
        
        # Now plot the 2d histogram
        img = ax.imshow(H,  interpolation='nearest', origin='lower',extent=[yedges[0], yedges[-1],xedges[0], xedges[-1]])
        colorbar = plt.colorbar(img, ax=ax)
        colorbar.set_label("Kcal / mol")
        
        # For normal histogram plot
        #plt.hist2d(np_arrays[0], np_arrays[1], bins=100) 
        #plt.colorbar()
                
        plt.ylim([-180,180])
        plt.xlim([-180,180])
        
        plt.title( title )
        plt.savefig(pp, format="pdf",dpi=150)
        pp.close()
        