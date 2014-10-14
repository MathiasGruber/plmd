#!/usr/bin/python
import os, shutil, sys
import tarfile
import caseAnalyses, globalAnalyses, plmd
import plmd.caseEmail

# Sending email modules
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

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
        
    # A wrapper for submitting the analysis to the HPC queue
    def submitForAnalysis( self, caseDir ):
        
        # User information
        self.printStage( "Creating submission script for server queue for folder: "+caseDir )
        
        # Python call on this dir
        sys.argv[0] = "plmd_analyse.py"
        sys.argv[2] = caseDir
        pythonCall = " ".join(sys.argv)   
        print ("PYTHON CALL:" + pythonCall)
        
        # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/analysis_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseDir  ). \
                              replace("[PYTHONCALL]", pythonCall )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseDir+"/submit_analysis.sh","w");        
        FILE.write( TEMP );
        FILE.close();

        # Submit the run
        os.system( "qsub "+caseDir+"/submit_analysis.sh" )
        
        # Remove submission file immidiately afterwards
        os.remove( caseDir+"/submit_analysis.sh" )
          
    # Main function handling analysis of a single case directory
    def analyseGlobal( self, searchDir ):
        
        # User information
        self.printStage( "Probing searchDir directory "+searchDir+" for data" )
    
        # Store data in this array
        dataArray = {}
        
        # Go through the directories
        for caseDir,dirs,files in os.walk( searchDir ):

            # Identify case structures
            if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs and "submit_run.sh" in files:
        
                # Only possible if it has a trajectory file
                if self.hasTrajectory( caseDir ):
                    
                    # Check if a data directory is present
                    if os.path.isdir( caseDir+"/analysis/data" ):
                        
                        # Go through the datafiles, add paths
                        for filename in os.listdir(  caseDir+"/analysis/data" ):
                            if filename not in dataArray:
                                dataArray[ filename ] = {'filepaths':[], 'caseLabels':[], 'filename':filename }
                            dataArray[ filename ]['filepaths'].append( caseDir+"/analysis/data/"+filename )
                            dataArray[ filename ]['caseLabels'].append( caseDir.split("/")[-1] )
        
        # Check that any data was found
        if len( dataArray ) > 0:
            
            # Create directory for the global analyses
            self.createFolder( "globalAnalyses" , True )            
            
            # Instantiate the handler for the analyses
            self.printStage( "Setting up analysis handler for: "+caseDir )
            handler = globalAnalyses.analysisHandler( self.configuration, dataArray )
            
            # Run all the analyses present in the handler
            handler.runAll()
            
            
        else:
            raise Exception("Not enough data was found to run global analysis")
        
                 
        
    # Main function handling analysis of a single case directory
    def analyseCase( self, caseDir ):
        
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
            self.createFolder( caseDir+"/analysis/structures" , True )
            self.createFolder( caseDir+"/analysis/structures/pca2" , True )
            self.createFolder( caseDir+"/analysis/structures/pca3" , True )
            self.createFolder( caseDir+"/analysis/structures/pca4" , True )
            
            # Instantiate the handler for the analyses
            self.printStage( "Setting up analysis handler for: "+caseDir )
            handler = caseAnalyses.analysisHandler( caseDir , self.configuration, self.num_files )
            
            # Run all ptraj analyses
            handler.runPtrajAnalysis()  
            
            # Run all the analyses present in the handler
            handler.runAll()
            
            # Create emailer
            emailObject = plmd.caseEmail.Setup( self.configuration )
            
            # Compress the analysis/plots folder
            folderToCompres = caseDir+"/analysis/plots"
            archieveName = caseDir+"/analysis/"+self.name+"-"+caseDir.split("/")[-1]
            emailObject.zipDirectory( archieveName , folderToCompres )
            
            # Email the compressed file to the user
            if self.noEmail == False:
                emailObject.emailFile( archieveName+".tar" , caseDir )
        
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
     