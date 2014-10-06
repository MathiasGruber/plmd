#!/usr/bin/python

try:
    
    # Standard libs
    import sys,os,traceback
    import argparse, ConfigParser
    
    # Disable use of display backend for matplotlib
    import matplotlib
    matplotlib.use('Agg')
    
    import plmd.caseAnalysis
    import getpass
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Scans through directories and analyses all production trajectories in case directories.
    """)
    
    # Only argument is the configFile
    parser.add_argument('configFile', help='Configuration file with details of the analysis')
    parser.add_argument('scanDir', help='Directory to scan for MD cases to submit to queue') 
    parser.add_argument('-noQueue', dest='noQueue',action='store_const', const="true", default="false", help='Do not run post-processing on the cluster queue') 
    parser.add_argument('-noMerge', dest='noMerge',action='store_const', const="true", default="false", help="Do not merge all found mdcrd files") 
    parser.add_argument('-noStrip', dest='noStrip',action='store_const', const="true", default="false", help="Do not strip water molecules in trajectory files")
    parser.add_argument('-noBlock', dest='noBlock',action='store_const', const="true", default="false", help="Do not run block averaging analysis")
    parser.add_argument('-noEnergy', dest='noEnergy',action='store_const', const="true", default="false", help="Do not run energy analysis")
    parser.add_argument('-noEmail', dest='noEmail',action='store_const', const="true", default="false", help="Do not attempt to send email")
    parser.add_argument('-emailPass', dest='emailPass',nargs="?", default="false", help="Email password to send results to")
    
    # Parse arguments
    args = parser.parse_args()  

    # Configuration
    config = ConfigParser.RawConfigParser()
    config.read( args.configFile )
    
    # Use command-line input to set config data
    config.set('analysisParameters', 'noMerge', args.noMerge )
    config.set('analysisParameters', 'noStrip', args.noStrip )
    config.set('analysisParameters', 'noBlock', args.noBlock )
    config.set('analysisParameters', 'noEnergy', args.noEnergy )
    config.set('analysisParameters', 'noEmail', args.noEmail )
    
    # Get email pass
    if args.noEmail == "false":
        if args.emailPass == "false":
            var = getpass.getpass("""
Using the smtp server specified in the configuration file,
this tool can email the results to your email account. 
Press 'n' to not use this feature\n
Please enter the password for your email, so as to use the smtp server.
The password is not permanently saved in any local or external files.
""")
            if var != 'n':
                config.set('emailConfiguration', 'emailPass', var )
                sys.argv.append( "-emailPass" )
                sys.argv.append( '"'+var+'"' )
        else:
            config.set('emailConfiguration', 'emailPass', args.emailPass )
        
    # The object handling analyses
    analyser = plmd.caseAnalysis.Analysis( config )    

    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):
        
        # Go through the directories
        for subdir,dirs,files in os.walk( args.scanDir ):

            # Identify case structures
            if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs and "submit_run.sh" in files:
                
                # Run analysis on this directory. 
                if args.noQueue == "true":
                    
                    # Directly in terminal
                    analyser.analyseDir( subdir )    
                else:
                    
                    # Submit to queue
                    sys.argv.append("-noQueue")
                    analyser.submitForAnalysis( subdir )
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