import os,sys
import MDAnalysis
import block, energy, bFactor, dihedral, timeCorr, endToEnd,CaToCaMap
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
        self.backbone = self.mdTrajectory.selectAtoms('protein and backbone')         
        self.simFrames = self.mdTrajectory.trajectory.numframes     
        self.simTime = self.mdTrajectory.trajectory.numframes * self.timestepSize * self.timestepPerFrame * 0.001
        
        # Print information about the trajectory for the user to see
        print "Number of frames in trajectory: "+str(self.mdTrajectory.trajectory.numframes)
        print "Number of timesteps per frame: "+str(self.timestepPerFrame)
        print "Timestep size: "+str(self.timestepSize)+"ps"
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
        
        # Do the dihedral angle specifications
        numOfResidues = self.backbone.numberOfResidues()
        dihedralTxt = ""
        if numOfResidues > 1:
            for i in range(1,numOfResidues):
                dihedralTxt += "\ndihedral phi_"+str(i)+" :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA :"+str(i+1)+"@C out "+caseDir+"/analysis/data/phi_"+str(i)
                dihedralTxt += "\ndihedral psi_"+str(i)+" :"+str(i)+"@N  :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N out "+caseDir+"/analysis/data/psi_"+str(i)
                dihedralTxt += "\ndihedral omega_"+str(i)+" :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA out "+caseDir+"/analysis/data/omega_"+str(i)
        
        # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/cpptraj_analysis.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseDir  ). \
                               replace("[DIHEDRALS]", dihedralTxt ). \
                               replace("[FIRSTRESI]", ":1" ). \
                               replace("[LASTRESI]", ":"+str(numOfResidues) )
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
        
    # Run a block analysis
    def dihedralAnalysis( self ):
        dihedral.runAnalysis( self.directory, self.backbone );
        
    # Time correlation
    def timeCorrelationAnalysis(self):
        timeCorr.runAnalysis( self.directory )
    
    # End to end distance
    def endtoendAnalysis(self):
        endToEnd.runAnalysis( self.directory )
    
    # Ca map
    def caMapAnalysis(self):
        CaToCaMap.runAnalysis(self.directory)