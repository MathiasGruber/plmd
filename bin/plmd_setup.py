#!/usr/bin/python

try:
    
    import argparse, ConfigParser, sys, traceback
    import plmd.caseSetup
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Sets up case directories and input files for Molecular Dynamics in Amber12.
    """)
    
    # Predefined ligands available in this package
    predefinedLigands = [ "PO4","HPO4","H2PO4","SO4","HSO4","H2SO4" ]
    
    # Only argument is the configFile
    parser.add_argument('configFile', help='a configuration file with details of the simulation')
    parser.add_argument('-pLigand', dest='ligand', help='Predefined Ligand from list: '+', '.join(predefinedLigands))
    parser.add_argument('-lCount', dest='ligandCount', help='Number of ligands to include.')
    parser.add_argument('-pPeptide', dest='peptide', help='3-letter notation string for automatic peptide generation with LEaP, e.g. "NSER GLY ALA GLY LYS CTHR"')
    parser.add_argument('-pCount', dest='peptideCount', help='Number of peptides to include. Only works for automatically LEaP generated peptides sequences.')
    parser.add_argument('--quiet', dest='quiet',action='store_const', const=True, default=False, help="Don't show status messages during setup.")
    
    # Parse arguments
    args = parser.parse_args()  
    
    # Check file can be opened
    f = open( args.configFile )
    
    # Open & parse the config file
    config = ConfigParser.RawConfigParser({'pLigand': None, 'pPeptide': None})
    config.read( args.configFile )
    
    # If a pLigand was specified, overwrite input file value
    if args.ligand:
        if args.ligand in predefinedLigands:
            print "Setting the ligand to: "+args.ligand
            config.set('inputFiles', 'pLigand', args.ligand )

        # If a count was specified, overwrite input file value
        if args.ligandCount:
            print "Setting the ligand count to: "+args.ligandCount
            config.set('inputFiles', 'ligandCount', args.ligandCount )        
        
    # If a pPeptide was specified, overwrite input file value
    if args.peptide:
        print "Setting the peptide to: "+args.peptide
        config.set('inputFiles', 'pPeptide', args.peptide )
        
        # If a count was specified, overwrite input file value
        if args.peptideCount:
            print "Setting the peptide count to: "+args.peptideCount
            config.set('inputFiles', 'peptideCount', args.peptideCount )
    
    # If quiet option was requested, then update the config var
    config.set('inputFiles', 'quiet', args.quiet )
    
    # Instantiate PLMD and pass the configuration file
    plmd = plmd.caseSetup.Setup( config )
    
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