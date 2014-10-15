import os, math, re
import MDAnalysis
import plmd
import energy, dihedral

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, config, dataFiles ):

        # Save config parameters
        self.load_config( config )

        # The directory we're working in 
        self.dataFiles = dataFiles

        # Load Trajectory of first case only. Just to get backbone structure
        self.mdTrajectory = MDAnalysis.Universe( self.dataFiles[ 'caseDirs' ][0]+"/md-files/peptide_nowat.prmtop", self.dataFiles[ 'caseDirs' ][0]+"/mergedResult.dcd") 
        self.backbone = self.mdTrajectory.selectAtoms('protein and backbone')    
        
    # Run all the analyses modules
    def runAll( self ):
    
        # Plot Energies
        if self.noEnergy == False:
            self.printStage( "Plotting energies")
            energy.runAnalysis( "Kinetic Energies", self.dataFiles[ "summary.EKTOT" ] );
            energy.runAnalysis( "Potential Energies", self.dataFiles[ "summary.EPTOT" ] );
            energy.runAnalysis( "Total Energies", self.dataFiles[ "summary.ETOT" ] );
            
        # Get dihedral angles
        dihedral.runAnalysis( self.dataFiles, self.backbone )          
        
     