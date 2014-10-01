#!/usr/bin/python
import argparse, sys
import plmd.configCreator
import traceback

# Argument parsing & help page
parser = argparse.ArgumentParser(description=
"""PLMD (Peptide Ligand Molecular Dynamics) 
Retrieves a default configuration file for use with PLMD
""")

# Start execution
try:
    # Parse arguments
    args = parser.parse_args()  
    

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