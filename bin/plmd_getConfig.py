#!/usr/bin/python

# Start execution
try:
    
    # Import needed modules
    import argparse, traceback, sys, ConfigParser, os
    import plmd.defaultConfig
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Retrieves a default configuration file for use with PLMD
    """)
    
    # Get a default config file
    config = plmd.defaultConfig.getDefaultConfig()
    dirName = os.getcwd().split(os.sep)[-1]
    config.set('plmd_settings', 'name', dirName )

    # Writing our configuration file to 'example.cfg'
    with open('defaultExample.cfg', 'wb') as configfile:
        config.write(configfile)
    
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