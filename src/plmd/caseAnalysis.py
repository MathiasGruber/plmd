#!/usr/bin/python
import os, shutil

class Analysis:

    def __init__(self, config):
        
        # Save config parameters
        self.email = config.get('analysisParameters', 'email')
        
        # Clear screen
        os.system("clear")        
        

        if self.quiet == False:
            var = raw_input("Please confirm this submission with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Configuration file was not confirmed')
            
    def analyseDir( self, caseDir ):
        print "Now going to analyse ",caseDir

    