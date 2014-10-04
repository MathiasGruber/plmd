#!/usr/bin/python
import os

class Analysis:

    def __init__(self, config):
        
        # Save config parameters
        self.email = config.get('analysisParameters', 'email')
        self.noMerge = config.get('analysisParameters', 'noMerge')
        self.noStrip = config.get('analysisParameters', 'noStrip')
        
        # Clear screen
        os.system("clear")        
            
    # Main function handling analysis of a single case directory
    def analyseDir( self, caseDir ):
        
        # User information
        self.printStage( "Now going to analyse: "+caseDir )
        
        # Merge Trajectories
        if self.noMerge == False:
            self.mergeTrajectories( caseDir )
        
        
    # A function for merging all the trajectories in a case fodler
    def mergeTrajectories( self, caseDir ):
        
        # Count the number of trajectory files
        path = caseDir+'/md-files/'
        num_files = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and ".mdcrd" in f and "equil" in f] ) 
        print("Number of trajectory files found: "+str(num_files))

        # Start creating ptraj script
        filename = caseDir+"/trajectoryMerge.ptraj";
        FILE = open(filename,"w");
        buffer = ""
        
        # Add all trajectory files to ptraj script
        for i in range(1,num_files):
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
        
        # Create output binpos file (~48% size of .mdcrd)
        buffer = buffer + """
        trajout """+caseDir+"""/mergedResult.binpos binpos nobox
        go
        """
        
        # Write & Close file
        FILE.write(buffer);
        FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+path+"/peptide.prmtop -i "+caseDir+"/trajectoryMerge.ptraj" )
        
    
    # A function for printing the current stage of the process to the user
    def printStage( self,info ):
        print "\n"+"="*len(info)
        print info
        print "="*len(info)+"\n"