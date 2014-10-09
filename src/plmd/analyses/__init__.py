import os, math
import MDAnalysis
import block, energy, bFactor, dihedral, timeCorr, endToEnd,CaToCaMap, RMSdMap, hbond
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
        
        # Set how many frames to skip in analysis
        # Ideally we'd want a maximum of 500 points/frames
        self.framesToSkip = 1
        if self.simFrames > 500:
            self.framesToSkip = math.floor( self.simFrames / 500. )
                
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
    def runPtrajAnalysis( self ):
        
        # Calculate factor for time-axis in ptraj analyses
        self.ptrajFactor = int(self.framesToSkip * self.timestepPerFrame * self.timestepSize)
        
        # Do the dihedral angle specifications
        numOfResidues = self.backbone.numberOfResidues()
        dihedralTxt = ""
        if numOfResidues > 1:
            for i in range(1,numOfResidues):
                dihedralTxt += "\ndihedral phi_"+str(i)+" :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA :"+str(i+1)+"@C out "+self.directory+"/analysis/data/phi_"+str(i)
                dihedralTxt += "\ndihedral psi_"+str(i)+" :"+str(i)+"@N  :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N out "+self.directory+"/analysis/data/psi_"+str(i)
                dihedralTxt += "\ndihedral omega_"+str(i)+" :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA out "+self.directory+"/analysis/data/omega_"+str(i)
        
        # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/cpptraj_analysis.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", self.directory  ). \
                               replace("[DIHEDRALS]", dihedralTxt ). \
                               replace("[FIRSTRESI]", ":1" ). \
                               replace("[LASTRESI]", ":"+str(numOfResidues) ). \
                               replace("[FRAMESKIP]", str(int(self.framesToSkip)) )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(self.directory+"/ccptraj_analysis.ptraj","w");        
        FILE.write( TEMP );
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+self.directory+"/md-files/peptide_nowat.prmtop -i "+self.directory+"/ccptraj_analysis.ptraj" )
        
    # Run all the analyses modules
    def runAll( self ):
    
        ## Run Analyses using MDAnalaysis module            
        ########################################
        
        # Block averaging analysis
        if self.noBlock == False:
            self.printStage( "Running block analysis for: "+self.directory )
            block.runAnalysis( self.directory, self.mdTrajectory, self.timestepSize );           
        
        ## Run analyses using cpptraj
        #############################
        
        # H-bond plotting
        hbond.runAnalysis( self.directory );            
        
        # Plot the B factor
        bFactor.runAnalysis( self.directory );
        
        # Plot angles
        dihedral.runAnalysis( self.directory, self.backbone , self.ptrajFactor);
        
        # Time correlation
        timeCorr.runAnalysis( self.directory , self.ptrajFactor)
        
        # Time correlation
        endToEnd.runAnalysis( self.directory , self.ptrajFactor)
        
        # C-alpha map
        CaToCaMap.runAnalysis(self.directory )
        
        # RMSd map
        RMSdMap.runAnalysis(self.directory, self.ptrajFactor)
        
        ## Run analyses using MD log data            
        #################################
            
        # Plot Energies
        if self.noEnergy == False:
            self.printStage( "Plotting energies: "+self.directory )
            energy.runAnalysis( self.directory );

     