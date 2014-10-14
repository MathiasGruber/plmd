import os, math
import MDAnalysis
import plmd
import energy

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, config, dataFiles ):

        # Save config parameters
        self.load_config( config )

        # The directory we're working in 
        self.dataFiles = dataFiles

        # Set how many frames to skip in analysis
        # Ideally we'd want a maximum of 500 points/frames
        #self.framesToSkip = 1
        #if self.simFrames > 500:
        #    self.framesToSkip = math.floor( self.simFrames / 500. )
        #self.ptrajFactor = int(self.framesToSkip * self.timestepPerFrame * self.timestepSize)
        
    # Run all the analyses modules
    def runAll( self ):
    
        # Plot Energies
        if self.noEnergy == False:
            self.printStage( "Plotting energies")
            energy.runAnalysis( "Kinetic Energies", self.dataFiles[ "summary.EKTOT" ] );
            energy.runAnalysis( "Potential Energies", self.dataFiles[ "summary.EPTOT" ] );
            energy.runAnalysis( "Total Energies", self.dataFiles[ "summary.ETOT" ] );

     