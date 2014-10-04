#!/usr/bin/python

# Start execution
try:
    
    import ConfigParser, argparse
    import traceback
    
    # Argument parsing & help page
    parser = argparse.ArgumentParser(description=
    """PLMD (Peptide Ligand Molecular Dynamics) 
    Retrieves a default configuration file for use with PLMD
    """)
    
    # The default config options
    config = ConfigParser.RawConfigParser()
    
    config.add_section('inputFiles')
    config.set('inputFiles', 'ligand', 'filenameLigand.mol2')
    config.set('inputFiles', 'ligandCount', '1')
    config.set('inputFiles', 'peptide', 'filenamePeptide.pdb')
    config.set('inputFiles', 'peptideCount', '1')
    config.set('inputFiles', 'cases', '10')
    
    config.add_section('simulationParameters')
    config.set('simulationParameters', 'forceField', 'leaprc.ff12SB')
    config.set('simulationParameters', 'qmCharge', '-2')
    config.set('simulationParameters', 'qmShake', '1')
    config.set('simulationParameters', 'qmTheory', 'PM6-D')
    config.set('simulationParameters', 'ntf', '2')
    config.set('simulationParameters', 'ntc', '2')
    config.set('simulationParameters', 'timestepSize', '0.002')
    config.set('simulationParameters', 'timestepNumber', '1000000')

    config.add_section('submissionParameters')
    config.set('submissionParameters', 'nodeControl', 'nodes=1:ppn=8')
    config.set('submissionParameters', 'wallClock', '72:00:00')
    config.set('submissionParameters', 'mdRuns', '50')

    config.add_section('analysisParameters')
    config.set('analysisParameters', 'email', 'nano.mathias@gmail.com')

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