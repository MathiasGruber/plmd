#!/usr/bin/env python

# Start execution
try:
    
    # Import needed modules
    import argparse, traceback, sys, os, re
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Goes through supplied cases and removes all MD runs above the specified number
    """)
    
    # Parameters
    parser.add_argument('scanDir', help='Directory to scan for MD cases to cut')
    parser.add_argument('cutoff', help='Delete all equilX.mdcrd and output files above this cutoff')
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Check that a directory was supplied
    if os.path.isdir( args.scanDir ):  
        
        # Check cutoff
        if args.cutoff > 1:        
        
            # Go through the directories
            for subdir,dirs,files in os.walk( args.scanDir ):
                
                # Check for the directories where we're deleting
                if "md-files" in subdir or "md-logs" in subdir or "pdb-files" in subdir:
                    for filename in files:
                        digits = re.findall(r'\d+', filename)
                        if digits:
                            if int(digits[0]) > int(args.cutoff):
                                toDelete = subdir+"/"+filename
                                print "Deleting: ", toDelete
                                os.system( "rm -rf "+toDelete )
           
        else:
            raise Exception("The specified cutoff is not valid")  
    else:
        raise Exception("The specified scanDir is not a directory")  
    
    
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