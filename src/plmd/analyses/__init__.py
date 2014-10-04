import MDAnalysis
import block, plmd

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, caseDir, config ):

        # Save config parameters
        self.load_config( config )

        # The directory we're working in 
        self.directory = caseDir

        # Load Trajectory with MDAnalysis
        self.mdTrajectory = MDAnalysis.Universe( caseDir+"/md-files/peptide_nowat.prmtop", caseDir+"/mergedResult.dcd")            
        self.simFrames = self.mdTrajectory.trajectory.numframes     
        self.simTime = self.mdTrajectory.trajectory.numframes * self.timestepSize
        
        # Print information about the trajectory for the user to see
        print "Number of frames in trajectory: "+str(self.mdTrajectory.trajectory.numframes)
        print "Total Simulation time: "+str( self.simTime )+"ns"

    # Run a block analysis
    def blockAnalysis( self ):
        block.runAnalysis( self.directory, self.mdTrajectory, self.timestepSize );