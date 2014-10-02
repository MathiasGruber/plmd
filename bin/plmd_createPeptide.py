#!/usr/bin/python

try:
    
    import argparse, sys
    import plmd.peptideConstructor as pepConstruct
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """
    Uses LEaP from AmberTools to create a pdb file with a linear peptide.
    """)
    
    # Only argument is the configFile
    parser.add_argument('i', help='3-letter notation string of peptide to generate with LEaP, e.g. "NSER GLY ALA GLY LYS CTHR"')
    parser.add_argument('o', help='where to output the pdb file')
    parser.add_argument('-n', dest='peptideCount', help='Number of peptides to include in .pdb-file (default: 1)', default=1)
    
    
    # Parse arguments
    args = parser.parse_args()  

    # Call the peptidce creator class    
    peptideCreator = pepConstruct.Creator()
    peptideCreator.createPeptide( args.i, int(args.peptideCount), args.o )
    
    
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