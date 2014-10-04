#!/usr/bin/python
import os
import analyses, plmd

# This is the overall Analysis class which merges trajectories,
# manages the analysis handler, and emails the final results
class Analysis (plmd.PLMD_module):

    def __init__(self, config):
        
        # Store the configuration data
        self.configuration = config        
        
        # Save config parameters
        self.load_config( config )  
        
        # Clear screen
        os.system("clear")        
            
    # Main function handling analysis of a single case directory
    def analyseDir( self, caseDir ):
        
        # User information
        self.printStage( "Analysis of case directory: "+caseDir )
        
        # Get the number of files
        self.num_files = self.getNumberOfFiles( caseDir+'/md-files/' ) 
        
        # Merge Trajectories
        if self.noMerge == False:
            self.mergeTrajectories( caseDir )
            
        # Run analyses if trajectory file is present
        if self.hasTrajectory( caseDir ):
            
            # Create the directory for all the postprocessing stuff
            self.createFolder( caseDir+"/analysis" , True )
            self.createFolder( caseDir+"/analysis/plots" , True )
            self.createFolder( caseDir+"/analysis/data" , True )
            
            # Instantiate the handler for the analyses
            self.printStage( "Setting up analysis handler for: "+caseDir )
            handler = analyses.analysisHandler( caseDir , self.configuration, self.num_files )
            
            # Run block analysis
            if self.noBlock == False:
                self.printStage( "Running block analysis for: "+caseDir )
                handler.blockAnalysis()
                
            # Plot Energies
            if self.noEnergy == False:
                self.printStage( "Plotting energies: "+caseDir )
                handler.energyAnalysis()
            
            # First compress the analysis/plots folder
            self.zipDirectory( caseDir+"/analysis/plots" )
            
            # Email the compressed file to the user
            self.emailFile( caseDir+"/analysis/plots.tar" )
        
    # A function for merging all the trajectories in a case fodler
    def mergeTrajectories( self, caseDir ):
        
        # Count the number of trajectory files
        path = caseDir+'/md-files/'
        print("Number of trajectory files found: "+str(self.num_files))

        # Start creating ptraj script
        filename = caseDir+"/trajectoryMerge.ptraj";
        FILE = open(filename,"w");
        buffer = ""
        
        # Add all trajectory files to ptraj script
        for i in range(1,self.num_files):
            buffer = buffer + """
            trajin """+path+"""equil""" + str(i)+ """.mdcrd"""
    
        # Center around first residue and in origin
        buffer = buffer + """
        center origin :1
        image origin center
        """
        
        # Strip water molecules
        if self.noStrip == False:
            buffer = buffer + "strip :WAT"
        
        # Create output binpos file 
        # Advantages of this format:
        # 1. ~48% size of .mdcrd
        # 2. ~same size as binpos, 
        # 3. Compatibility with MDAnalysis module
        buffer = buffer + """
        trajout """+caseDir+"""/mergedResult.dcd charmm nobox
        go
        """
        
        # Write & Close file
        FILE.write(buffer);
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+path+"/peptide.prmtop -i "+caseDir+"/trajectoryMerge.ptraj" )
        
    # Get number of files
    def getNumberOfFiles( self, path ):        
        return len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and ".mdcrd" in f and "equil" in f] ) 
        
    # A function for checking whether a binpos file exists
    def hasTrajectory( self, caseDir ):
        if os.path.isfile( caseDir+"/mergedResult.dcd" ):
            return True
            
    # Zip a given directory
    def zipDirectory( self, caseDir ):
        self.printStage("Compressing: "+caseDir)
        
    # Email a given file to the user
    def emailFile( self, filepath ):
        self.printStage("Emailing the file: "+filepath)
    