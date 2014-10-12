#!/usr/bin/python
import os
import plmd 

class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Load the config file
        self.load_config( config )        
        
        # Confirm with user
        self.printStage("Step 1: Starting up PLMD. Submission file parameters:")

        print "\n== Submission Parameters"
        print "========================"
        print "submissionName: " + self.name
        print "nodeControl: " + self.nodeControl
        print "wallClock: " + self.wallClock
        print "mdRuns: " + self.mdRuns

        # Confirmation from user
        if self.quiet == False:
            var = raw_input("\nPlease confirm these submission parameters with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Submission was not confirmed')

    
    # Create submission file for submitting case to HPC queue
    def hpcCreateSubmission( self, caseName ):
        
        # User information
        self.printStage( "Stage 2, Case: "+caseName+". Creating HPC submission files" )          
        caseID = caseName.split("/")[-1] 
        
         # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/explicit_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseName  ). \
                              replace("[NAME]", self.name+"_"+caseID  ). \
                              replace("[CPUCONTROL]", self.nodeControl ). \
                              replace("[WALLCLOCK]", self.wallClock ). \
                              replace("[MDRUNS]", self.mdRuns )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseName+"/submit_run.sh","w");        
        FILE.write( TEMP );
        FILE.close();
        
        print "Create submission file: "+caseName+"/submit_run.sh"
    # Submit to HPC cluster
    def hpcSubmission( self, caseName ):
        
        # User information
        self.printStage( "Stage 3, Case: "+caseName+". Submitting to HPC" )          
        
        # Do submission        
        os.system( "qsub "+caseName+"/submit_run.sh" )
             
            