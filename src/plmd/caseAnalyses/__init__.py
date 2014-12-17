import os, math, sys
import MDAnalysis
import block, energy, bFactor, dihedral, timeCorr, endToEnd,CaToCaMap, RMSdMap, hbond, pca, hierCluster, dbscanCluster, RMSdFrequency, mmpbsa
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
        
        # If aMD files are found, merge them into a weights file
        if os.path.isfile(caseDir+"/md-logs/aMD1.log"): 
            self.printStage( "Creating aMD weights file" )
            aMDfiles = len([f for f in os.listdir(caseDir+"/md-logs/") if os.path.isfile(os.path.join(caseDir+"/md-logs/", f)) and "aMD" in f and "out" not in f and ".log" in f] )
            print "Found "+str(aMDfiles)+" aMD log files" 
            
            # Get the files as per specified at: http://mccammon.ucsd.edu/computing/amdReweighting/
            for i in range(1,aMDfiles+1):
                os.system("awk 'NR>3' "+caseDir+"/md-logs/aMD"+str(i)+".log | awk '{print ($8+$7)/(0.001987*300)\" \" $2 \" \" ($8+$7)}' > "+caseDir+"/md-logs/weights"+str(i)+".dat")
                
            # Merge all of them
            buffer = ""
            stepCounter = 0
            for i in range(1,aMDfiles+1):
                fileCounter = 0
                with open(caseDir+"/md-logs/weights"+str(i)+".dat", "r") as fi:
                    for line in fi:
                        temp = line.split()
                        increase = int(temp[1])-fileCounter
                        stepCounter += increase
                        fileCounter = int(temp[1])
                        buffer += "\n"+temp[0]+" "+str(stepCounter)+" "+temp[2]
                        
            # Write the weights time
            with open(caseDir+"/md-logs/weights.dat", "w") as fo:
                fo.write( buffer )
            
        # Run perl script from AMBER to get data from log files
        buffer = ""
        for i in range(1,num_files):
            if self.config.amdEnabled == True:
                if os.path.isfile(caseDir+"/md-logs/outAMD"+str(i)+".log"): 
                    buffer = buffer + " "+caseDir+"/md-logs/outAMD"+str(i)+".log"
                elif os.path.isfile(caseDir+"/md-logs/outMD"+str(i)+".log"):
                    buffer = buffer + " "+caseDir+"/md-logs/outMD"+str(i)+".log"    
            else:
                if os.path.isfile(caseDir+"/md-logs/outMD"+str(i)+".log"):
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
        
        # We run ptraj twice, once going through all frames (for analyses requiring that)
        # and once limited the trajectory to just enough points for a good plot
        for ptrajType in ["full", "short"]:        
        
            # Calculate factor for time-axis in ptraj analyses
            self.ptrajFactor = int(self.framesToSkip * self.config.timestepPerFrame * self.config.timestepSize)
            
            # Create new submission file
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/cpptraj_analysis_"+ptrajType+".txt", 'r')
            TEMP = TEMPLATE.read().replace("[FOLDER]", self.directory  ). \
                                   replace("[DIHEDRALS]", dihedralTxt ). \
                                   replace("[FIRSTRESI]", "1" ). \
                                   replace("[SIEVE]", str(self.sieveValue) ). \
                                   replace("[MINPOINTS]", str(self.config.dbscanMinPoints) ). \
                                   replace("[DBSCANEPS]", str(self.config.dbscanEps) ). \
                                   replace("[LASTRESI]", str(numOfResidues) ). \
                                   replace("[LASTID]", str( numOfResidues + self.config.ligandCount ) ). \
                                   replace("[FRAMESKIP]", str(int(self.framesToSkip)) )
            TEMPLATE.close()
                                  
            # Move k-dist from dbScan cluster to data directory
            os.system("mv Kdist.*.dat "+self.directory+"/analysis/data/")                                  
                                  
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
        try:
            if self.config.noBlock == False:
                self.printStage( "Running block analysis for: "+self.directory )
                block.runAnalysis( self.directory, self.mdTrajectory, self.timePerFrame );   
        except Exception as e:
            print "Failed block analysis",e
        
        ## Run analyses using cpptraj
        #############################

         # MMPBSA Energies
        try:
            mmpbsa.runAnalysis( self.directory, self.timePerFrame  )
        except Exception as e:
            print "Failed mmpbsa analysis",e

        # RMSd Frequencies
        try:
            RMSdFrequency.runAnalysis( self.directory )
        except Exception as e:
            print "Failed block analysis",e
        
        # Hieragglo cluster
        try:
            hierCluster.runAnalysis( self.directory )
        except Exception as e:
            print "Failed hier cluster analysis",e
        
        # dbscan cluster
        try:
            dbscanCluster.runAnalysis( self.directory, self.config.dbscanEps, self.config.dbscanMinPoints )
        except Exception as e:
            print "Failed dbscan cluser analysis",e
        
        # H-bond plotting
        try:
            hbond.runAnalysis( self.directory );            
        except Exception as e:
            print "Failed hbond analysis",e
        
        # RMSd map
        try:
            RMSdMap.runAnalysis(self.directory, self.ptrajFactor)        
        except Exception as e:
            print "Failed rmsdmap analysis",e
        
        # Plot angles
        try:
            dihedral.runAnalysis( self.directory, self.backbone , self.timePerFrame );
        except Exception as e:
            print "Failed dihedral analysis",e
              
        # Run the PCA analysis
        try:
            pca.runAnalysis( self.directory , self.mdTrajectory )        
        except Exception as e:
            print "Failed pca analysis",e
        
        # Plot the B factor
        try:
            bFactor.runAnalysis( self.directory );
        except Exception as e:
            print "Failed bfactor analysis",e
        
        # Time correlation
        try:
            timeCorr.runAnalysis( self.directory , self.ptrajFactor)
        except Exception as e:
            print "Failed timecorr analysis",e
        
        # Time correlation
        try:
            endToEnd.runAnalysis( self.directory , self.ptrajFactor)
        except Exception as e:
            print "Failed endtoend analysis",e
        
        # C-alpha map
        try:
            CaToCaMap.runAnalysis(self.directory )
        except Exception as e:
            print "Failed CaToCa map analysis",e
        
        
        ## Run analyses using MD log data            
        #################################
            
        # Plot Energies
        try:
            if self.config.noEnergy == False:
                self.printStage( "Plotting energies: "+self.directory )
                energy.runAnalysis( self.directory );
        except Exception as e:
            print "Failed energy analysis",e
     