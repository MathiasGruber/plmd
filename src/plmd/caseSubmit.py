#!/usr/bin/python
import os
import plmd 

class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Load the config file
        self.config = plmd.PLMD_Config( config ) 
        
        # Confirm with user
        self.printStage("Step 1: Starting up PLMD. Submission file parameters:")

        # Add GPU to nodecontrol if applicable
        if self.config.gpuEnabled == True:
            self.config.nodeControl += ":gpus="+str(self.config.gpuCores)

        print "\n== Submission Parameters"
        print "========================"
        print "submissionName: " + self.config.name
        print "nodeControl: " + self.config.nodeControl
        print "wallClock: " + self.config.wallClock
        print "mdRuns: " + self.config.mdRuns

        # Confirmation from user
        if self.config.quiet == False:
            var = raw_input("\nPlease confirm these submission parameters with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Submission was not confirmed')

    
    # Create submission file for submitting case to HPC queue
    def hpcCreateSubmission( self, caseName ):
        
        # User information
        self.printStage( "Stage 2, Case: "+caseName+". Creating HPC submission files" )          
        caseID = caseName.split("/")[-1] 
        
         # Create new submission file
        TEMPLATE = ""
        if self.config.gpuEnabled == True:
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/explicit_gpu_submit.txt", 'r')
        else:
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/explicit_submit.txt", 'r')
            
        # Replace stuff within
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseName  ). \
                              replace("[NAME]", self.config.name+"_"+caseID  ). \
                              replace("[CPUCONTROL]", self.config.nodeControl ). \
                              replace("[WALLCLOCK]", self.config.wallClock ). \
                              replace("[MDRUNS]", self.config.mdRuns )
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
             
            