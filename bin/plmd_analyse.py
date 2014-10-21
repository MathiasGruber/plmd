#!/usr/bin/env python

try:
    
    # Standard libs
    import sys,os,traceback
    import argparse, ConfigParser
    
    # Require python 2.7 or higher
    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
        print("This script requires Python version 2.7")
        sys.exit(1)    
    
    # Disable use of display backend for matplotlib
    import matplotlib
    matplotlib.use('Agg')
    matplotlibVersion = matplotlib.__version__.split(".")
    if int(matplotlibVersion[0]) < 0 or int(matplotlibVersion[1]) < 4:
        print "This script requires matplotlib 1.4.0"
        print "MatPlotLIb version: ",matplotlib.__version__, matplotlib.__file__   
        sys.exit(1)    
    
    # Set plotting configuration
    import matplotlib.pyplot as plt

    plt.rcParams['axes.color_cycle'] = ["348ABD", "7A68A6", "A60628", "467821", "CF4457", "188487", "E24A33","348ABD", "7A68A6", "A60628", "467821", "CF4457", "188487", "E24A33"]
    
    
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
    parser.add_argument('-noFTP', dest='noFTP',action='store_const', const="true", default="false", help="Do not attempt to send results to FTP")
    parser.add_argument('-ftpPass', dest='ftpPass',nargs="?", default="false", help="FTP password to send results to")
    parser.add_argument('-doGlobal', dest='doGlobal',action='store_const', const="true", default="false", help='Perform global analyses of selected cases. This assumes that all local analyses have already been performed, and uses that data') 
    
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
    config.set('analysisParameters', 'noFTP', args.noFTP )
    
    # Get ftp pass
    if args.noFTP == "false":
        if args.ftpPass == "false":
            var = getpass.getpass("""
Using the ftp server specified in the configuration file,
this tool can send the results to that ftp account. 
Press 'n' to not use this feature\n
Please enter the password for your account, so as to use the ftp server.
The password is not permanently saved in any local or external files.
""")
            if var != 'n':
                config.set('ftpConfiguration', 'ftpPass', var )
                sys.argv.append( "-ftpPass" )
                sys.argv.append( '"'+var+'"' )
        else:
            config.set('ftpConfiguration', 'ftpPass', args.ftpPass )
        
    # The object handling analyses
    analyser = plmd.caseAnalysis.Analysis( config )    


    
    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):
        
        # Check if global
        if args.doGlobal == "true":
            
            # Run analysis on this directory. 
            if args.noQueue == "true":
                analyser.analyseGlobal( args.scanDir )    
            else:
                sys.argv.append("-noQueue")
                analyser.submitForAnalysis( args.scanDir )
            
        else:        
        
            # Go through the directories
            for subdir,dirs,files in os.walk( args.scanDir ):
    
                # Identify case structures
                if "in_files" in dirs and "md-files" in dirs and "md-logs" in dirs and "pdb-files" in dirs and "submit_run.sh" in files:
                    
                    # Run in terminal, or submit to server queue
                    if args.noQueue == "true":
                        analyser.analyseCase( subdir )    
                    else:
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