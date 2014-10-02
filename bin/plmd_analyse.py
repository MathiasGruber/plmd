#!/usr/bin/python

try:
    
    import argparse, ConfigParser, sys, os
    import plmd.baseClass
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    is designed to aid in setting up and running MD simulations of small
    peptides and ionic ligands using Amber and AmberTools MD packages. 
    """)
    
    # Only argument is the configFile
    parser.add_argument('configFile', help='a configuration file with details of the simulation')    
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Check file can be opened
    f = open( args.configFile )
    
    # Open & parse the config file
    config = ConfigParser.RawConfigParser()
    config.read( args.configFile )
    
    # Instantiate PLMD and pass the configuration file
    plmd = plmd.baseClass.PLMD_base( config )
    
    # Setup all the cases
    plmd.setupCases()
    
    
except IOError as e:
    print 'Cannot open file: ',e
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