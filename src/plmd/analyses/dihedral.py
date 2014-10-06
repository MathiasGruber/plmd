#!/usr/bin/python
import MDAnalysis
import plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , backbone ):

    # Retrieve info on the peptide
    resNames = backbone.resnames()
    
    # Go through each residue connection
    for i in range( 1, len(resNames) ):
        
        # User info
        print "Plotting dihedrals for residue: "+str(i)
        
        # Plot command
        myPlot.plotData( 
        caseDir , 
        "Dihedral Angles. Res ID: "+str(i)+", Residue name: "+str(resNames[i]), 
        ["$\Psi_"+str(i)+"$","$\Phi_"+str(i)+"$"],
        ["psi_"+str(i),"phi_"+str(i)] , 
        "Angle (Degrees)")
 