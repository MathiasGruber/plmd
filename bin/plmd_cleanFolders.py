#!/usr/bin/env python

try:

    # Import modules
    import argparse, sys, os, re
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Clears up unimportant logging data from MD runs, but retains the initial case structure & results.
    """)    
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Unimportant files
    dumbFiles = ["^leap.log$", "^.+\.e\d+$", "^.+\.o\d+$", "^Error.+\.err$", "^Output.+\.out$", "^.+\.tar$", "^LEaP_setup.log$", "^ccptraj_analysis_.+\.ptraj", "^MMPBSATEMP.+", "reweight_variable.log", "mdinfo"]    
    
    # Do a walkthrough
    for subdir,dirs,files in os.walk( "." ):
        for filename in files:
            for check in dumbFiles:
                if re.search(check, filename) != None:
                    os.remove( subdir+"/"+filename )    
                    print "Deleting file: "+subdir+"/"+filename
    
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