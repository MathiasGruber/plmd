import os, math
import MDAnalysis
import block, energy, bFactor, dihedral, timeCorr, endToEnd,CaToCaMap, RMSdMap, hbond, pca, hierCluster, dbscanCluster, RMSdFrequency
import plmd

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, caseDir, config, num_files ):

        # Save config parameters
        self.config = config

        # The directory we're working in 
        self.directory = caseDir

        # Load Trajectory with MDAnalysis
        self.mdTrajectory = MDAnalysis.Universe( caseDir+"/md-files/peptide_nowat.prmtop", caseDir+"/mergedResult.dcd") 
        self.backbone = self.mdTrajectory.selectAtoms('protein and backbone')         
        self.simFrames = self.mdTrajectory.trajectory.numframes     
        self.simTime = self.mdTrajectory.trajectory.numframes * self.config.timestepSize * self.config.timestepPerFrame * 0.001
        
        self.timePerFrame = self.config.timestepPerFrame * self.config.timestepSize        
        
        # Set how many frames to skip in analysis
        # Ideally we'd want a maximum of 500 points/frames
        self.framesToSkip = 1
        if self.simFrames > 500:
            self.framesToSkip = math.floor( self.simFrames / 500. )
                
        # Print information about the trajectory for the user to see
        print "Number of frames in trajectory: "+str(self.mdTrajectory.trajectory.numframes)
        print "Number of timesteps per frame: "+str(self.config.timestepPerFrame)
        print "Timestep size: "+str(self.config.timestepSize)+"ps"
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

        # Do the dihedral angle specifications
        numOfResidues = self.backbone.numberOfResidues()
        dihedralTxt = ""
        if numOfResidues > 1:
            for i in range(1,numOfResidues):
                dihedralTxt += "\ndihedral phi_"+str(i)+" :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA :"+str(i+1)+"@C out "+self.directory+"/analysis/data/phi_"+str(i)
                dihedralTxt += "\ndihedral psi_"+str(i)+" :"+str(i)+"@N  :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N out "+self.directory+"/analysis/data/psi_"+str(i)
                dihedralTxt += "\ndihedral omega_"+str(i)+" :"+str(i)+"@CA :"+str(i)+"@C  :"+str(i+1)+"@N  :"+str(i+1)+"@CA out "+self.directory+"/analysis/data/omega_"+str(i)
        
        # Sieve limits. Aim for 10000 frames in cluster analysis
        self.sieveValue = math.ceil( self.simFrames / self.config.clusterFrames )
        self.dbscanMinpoints = math.ceil(self.config.clusterFrames / 1000)
        if self.dbscanMinpoints < 10:
            self.dbscanMinpoints = 10
        
        
        # We run ptraj twice, once going through all frames (for analyses requiring that)
        # and once limited the trajectory to just enough points for a good plot
        for ptrajType in ["full", "short", "global_pca", "global_rmsd"]:        
        
            # Calculate factor for time-axis in ptraj analyses
            self.ptrajFactor = int(self.framesToSkip * self.config.timestepPerFrame * self.config.timestepSize)
            
            # Create new submission file
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/cpptraj_analysis_"+ptrajType+".txt", 'r')
            TEMP = TEMPLATE.read().replace("[FOLDER]", self.directory  ). \
                                   replace("[DIHEDRALS]", dihedralTxt ). \
                                   replace("[FIRSTRESI]", "1" ). \
                                   replace("[SIEVE]", str(self.sieveValue) ). \
                                   replace("[MINPOINTS]", str(self.dbscanMinpoints) ). \
                                   replace("[LASTRESI]", str(numOfResidues) ). \
                                   replace("[LASTID]", str( numOfResidues + self.config.ligandCount ) ). \
                                   replace("[FRAMESKIP]", str(int(self.framesToSkip)) )
            TEMPLATE.close()
                                  
            # Write the submission file
            FILE = open(self.directory+"/ccptraj_analysis_"+ptrajType+".ptraj","w");        
            FILE.write( TEMP );
            FILE.close();
            
            # Run the cpptraj utility
            os.system( "$AMBERHOME/bin/cpptraj -p "+self.directory+"/md-files/peptide_nowat.prmtop -i "+self.directory+"/ccptraj_analysis_"+ptrajType+".ptraj" )
        
    # Run all the analyses modules
    def runAll( self ):
    
        ## Run Analyses using MDAnalaysis module            
        ########################################
        
        # Block averaging analysis
        if self.config.noBlock == False:
            self.printStage( "Running block analysis for: "+self.directory )
            block.runAnalysis( self.directory, self.mdTrajectory, self.timePerFrame );           
        
        ## Run analyses using cpptraj
        #############################

        # RMSd Frequencies
        RMSdFrequency.runAnalysis( self.directory )

        # Hieragglo cluster
        hierCluster.runAnalysis( self.directory )
        
        # dbscan cluster
        dbscanCluster.runAnalysis( self.directory )

        # H-bond plotting
        hbond.runAnalysis( self.directory );            
        
        # RMSd map
        RMSdMap.runAnalysis(self.directory, self.ptrajFactor)        
        
        # Plot angles
        dihedral.runAnalysis( self.directory, self.backbone , self.timePerFrame );
                
        # Run the PCA analysis
        pca.runAnalysis( self.directory , self.mdTrajectory )        
        
        # Plot the B factor
        bFactor.runAnalysis( self.directory );
        
        # Time correlation
        timeCorr.runAnalysis( self.directory , self.ptrajFactor)
        
        # Time correlation
        endToEnd.runAnalysis( self.directory , self.ptrajFactor)
        
        # C-alpha map
        CaToCaMap.runAnalysis(self.directory )
        
        
        ## Run analyses using MD log data            
        #################################
            
        # Plot Energies
        if self.config.noEnergy == False:
            self.printStage( "Plotting energies: "+self.directory )
            energy.runAnalysis( self.directory );

     