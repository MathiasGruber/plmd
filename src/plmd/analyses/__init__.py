import os
import MDAnalysis
import block, energy, bFactor
import plmd

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, caseDir, config, num_files ):

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
        
        # Run perl script from AMBER to get data from log files
        buffer = ""
        for i in range(1,num_files):
            buffer = buffer + " "+caseDir+"/md-logs/outMD"+str(i)+".log"
        print buffer
        os.system("perl $PLMDHOME/src/perl/process_mdout.perl "+buffer)

        # Move all summary files to analysis/data folder
        os.system("mv summary* "+caseDir+"/analysis/data")

    # Create and run ptraj file
    def runPtrajAnalysis( self , caseDir ):
        
        # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/cpptraj_analysis.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseDir  )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseDir+"/ccptraj_analysis.ptraj","w");        
        FILE.write( TEMP );
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+caseDir+"/md-files/peptide_nowat.prmtop -i "+caseDir+"/ccptraj_analysis.ptraj" )
        

    # Run a block analysis
    def blockAnalysis( self ):
        block.runAnalysis( self.directory, self.mdTrajectory, self.timestepSize );
        
    # Run a block analysis
    def energyAnalysis( self ):
        energy.runAnalysis( self.directory );
        
        # Run a block analysis
    def bFactorAnalysis( self ):
        bFactor.runAnalysis( self.directory );