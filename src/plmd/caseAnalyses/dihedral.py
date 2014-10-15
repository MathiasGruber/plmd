#!/usr/bin/python
import MDAnalysis
import plmd.plotData as myPlot

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
        
        # Print time corrected files for later plotting
        for timeCorrection in [ psiPath, phiPath ]:
            with open( timeCorrection , "r" ) as fl:
                lines = fl.readlines()
            with open( timeCorrection+"_timeCorrected" , "w") as fo:
                n = 0
                for line in lines:
                    if n > 0:
                        temp = line.split()
                        fo.write( str(float(temp[0])*timeFactor) + "\t" + temp[1] +"\n")
                    else:
                        fo.write( line )
                    n += 1
        
        # Plot command
        myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "Dihedral Angles. Res ID: "+str(i)+", Residue name: "+str(resNames[i]), 
        ["$\Psi_"+str(i)+"$","$\Phi_"+str(i)+"$"],
        [ psiPath , phiPath ] , 
        "Angle (Degrees)", 
        xFactor = timeFactor )
 