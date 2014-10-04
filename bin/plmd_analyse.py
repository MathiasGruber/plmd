#!/usr/bin/python

try:
    
    import argparse, ConfigParser, sys, os
    import plmd.caseAnalysis
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Scans through directories and analyses all production trajectories in case directories.
    """)
    
    # Only argument is the configFile
    parser.add_argument('configFile', help='Configuration file with details of the analysis')
    parser.add_argument('scanDir', help='Directory to scan for MD cases to submit to queue') 
    parser.add_argument('-noMerge', dest='noMerge',action='store_const', const=True, default=False, help="Do not merge all found mdcrd files") 
    parser.add_argument('-noStrip', dest='noStrip',action='store_const', const=True, default=False, help="Do not strip water molecules in trajectory files") 
    
    # Parse arguments
    args = parser.parse_args()  

    # Configuration
    config = ConfigParser.RawConfigParser()
    config.read( args.configFile )
    
    # Use command-line input to set config data
    config.set('analysisParameters', 'noMerge', args.noMerge )
    config.set('analysisParameters', 'noStrip', args.noStrip )
    
    # The object handling analyses
    analyser = plmd.caseAnalysis.Analysis( config )    

    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):
        
        # Go through the directories
        for subdir,dirs,files in os.walk( args.scanDir ):

            # Identify case structures
            if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs and "submit_run.sh" in files:
                
                # Submit to HPC cluster
                print "Identified case structure for directory: "+subdir+". Submitting for analysis."
                analyser.analyseDir( subdir )    
    
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