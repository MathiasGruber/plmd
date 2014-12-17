#!/usr/bin/env python

try:
    
    import argparse, os, sys
    import traceback, ConfigParser
    import plmd.caseSubmit
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Scans supplied directory for PLMD cases and calculates the binding energy using MMPBSA.py
    """)
    
    # Only argument is the configFile
    parser.add_argument('configFile', help='a configuration file with details of the simulation')
    parser.add_argument('scanDir', help='Directory to scan for MD cases to submit to queue')    
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):
        
        # Check file can be opened
        f = open( args.configFile )

        # Get configuration
        config = ConfigParser.RawConfigParser()
        config.read( args.configFile )        
        
        # Instantiate PLMD and pass the configuration file
        plmd = plmd.caseSubmit.Setup( config )
        
        # Go through the directories
        for subdir,dirs,files in os.walk( args.scanDir ):

            # Identify case structures
            if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs:
                
                # Create submission file
                plmd.hpcMMPBSASubmissionCreate( subdir )
                
                # Create submission file
                plmd.hpcMMPBSASubmission( subdir )

                
    else:
        raise Exception("The specified scanDir is not a directory")           
                
    
except ImportError as e:
    print e
except ImportError as e:
    print e
except Exception as e:
    
    # Get into
    exc_type, exc_value, exc_traceback = sys.exc_info()
    
    # Show traceback
    print "\n\n====================="
    print "== ERROR TRACEBACK =="
    print "====================="
    for err in traceback.format_exception(exc_type, exc_value,exc_traceback):
        print err
    print "Script execution was stopped."