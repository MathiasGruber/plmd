#!/usr/bin/env python

try:
    
    import argparse, os, sys
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Scans supplied directory for PLMD cases and submits all found cases using their submission script
    """)
    
    # Only argument is the configFile
    parser.add_argument('scanDir', help='Directory to scan for MD cases to submit to queue')    
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):
        
        # Go through the directories
        for subdir,dirs,files in os.walk( args.scanDir ):

            # Identify case structures
            if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs and "submit_run.sh" in files:
                
                # Submit to HPC cluster
                print "Identified case structure for directory: "+subdir+". Now submitting to HPC queue."
                os.system( "qsub "+subdir+"/submit_run.sh" )
    
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