#!/usr/bin/python
import argparse, sys
import traceback

# Argument parsing & help page
parser = argparse.ArgumentParser(description=
"""PLMD (Peptide Ligand Molecular Dynamics) 
Clears up unimportant logging data from MD runs, but retains the initial case structure & results.
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