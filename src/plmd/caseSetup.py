#!/usr/bin/python
import os, shutil
from copy import deepcopy
import structureManipulation, peptideConstructor, plmd 

class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Load the config file
        self.load_config( config )        
        
        # Clear screen
        os.system("clear")        
        
        # Ensure integrity of input files
        self.setupInputFiles()
        
        # Setup quick links to tools in amber
        self.xLEaP = self.AMBERHOME+"/bin/tleap -f "

        # Confirm with user
        self.printStage("Step 1: Starting up PLMD. Configuration file:")

        print "\n== Input Files"
        print "=============="
        print "ligand: " + self.ligand
        print "ligand count: " + str(self.ligandCount)
        print "peptide: " + self.peptide
        print "peptide count: " + str(self.peptideCount)
        print "cases: " + str(self.cases)

        print "\n== Simulation Parameters"
        print "========================"
        print "Forcefield: " + self.ff
        print "qmCharge: " + self.qmCharge 
        print "qmShake: " + self.qmShake
        print "qmTheory: " + self.qmTheory
        print "ntf: " + self.ntf 
        print "ntc: " + self.ntc
        print "timestepSize: " + str(self.timestepSize)
        print "timestepNumber: " + self.timestepNumber

        print "\n== Submission Parameters"
        print "========================"
        print "nodeControl: " + self.nodeControl
        print "wallClock: " + self.wallClock
        print "mdRuns: " + self.mdRuns

        if self.quiet == False:
            var = raw_input("Please confirm this submission with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Configuration file was not confirmed')

    # In case a special peptide was specified, or a predifned ligand was chosen
    def setupInputFiles(self):
    
        # User information
        self.printStage( "Stage 0. Ensuring integrity of input files" )
        
        # Create directory
        self.createFolder ( "predefinedInput" )       
               
        # Check if a peptide should be made
        if self.pPeptide != None:

            # User info
            print "Now creating peptide with LEaP"
            
            # Create & save the peptide
            peptideCreator = peptideConstructor.Creator()
            peptideCreator.createPeptide( self.pPeptide, self.peptideCount, "predefinedInput/peptides.pdb" )
            
            # Overwrite default configs
            self.peptide = "predefinedInput/peptides.pdb"
              
        # Check if we should use predefined ion from package
        if self.pLigand != None:
            
            # User info
            print "Now importing the predefined ion to be user to library."
            
            # Check that the file exists in pre-defined library
            if os.path.isfile( self.PLMDHOME+"/src/ligands/"+self.pLigand+".mol2" ):
                
                # Copy and set
                shutil.copy( self.PLMDHOME+"/src/ligands/"+self.pLigand+".mol2", "predefinedInput/"+self.pLigand+".mol2" )
                self.ligand = "predefinedInput/"+self.pLigand+".mol2"
                
            else:
                
                # Raise error
                raise Exception("Tried to include invalid pre-defined ligand.")                

    # Main functon for setting up all the cases
    def setupCases(self):
        
        # Get the coordinates & resnames of peptide & ligand. Passed by reference.
        structureManipulation.processStructureFile( self.peptide , "pdb", self.peptideCoordinates, self.peptideResnames )
        structureManipulation.processStructureFile( self.ligand , "mol2", self.ligandCoordinates, self.ligandResnames )
        
        # Get center of masses
        self.peptideCenterOfMass = structureManipulation.centerOfMass( self.peptideCoordinates )
        self.ligandCenterOfMass = structureManipulation.centerOfMass( self.ligandCoordinates )
        
        # Create main folder for project
        self.createMainFolder()        
        
        # Go through all the cases to setup
        for i in range(0,self.cases):
            
            # Create folder structure
            self.creatCaseFolder( str(i) )            
            
            # User information
            self.printStage( "Stage 4, Case: "+str(i)+". Calculation ion translation vector" )
            
            # Calculate a random translation for the ligand (pass by value)
            # LEaP performs the translation using this vector
            self.ligandsForLEaP[:] = []
            for x in range(0, self.ligandCount):
                self.ligandsForLEaP.append( 
                    structureManipulation.calcIonPosition( 
                        deepcopy(self.peptideCenterOfMass), 
                        deepcopy(self.peptideCoordinates),
                        deepcopy(self.ligandCenterOfMass), 
                        deepcopy(self.ligandCoordinates),
                        self.ligandsForLEaP
                    )
                )
            
            # Create LEaP input File for case
            self.leapCreateInput( str(i) )
            
            # Run the LEaP file
            self.leapRunInput( str(i) )
            
            # Create input files for amber run
            self.amberCreateInput( str(i) )
            
            # Create files for HPC submission
            self.hpcCreateSubmission( str(i) )
    
    # Create submission file for submitting case to HPC queue
    def hpcCreateSubmission( self, caseName ):
        
        # User information
        self.printStage( "Stage 7, Case: "+caseName+". Creating HPC submission files" )          
        
         # Create new submission file
        TEMPLATE = open( self.PLMDHOME+"/src/templates/explicit_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", "cases/"+caseName  ). \
                              replace("[CPUCONTROL]", self.nodeControl ). \
                              replace("[WALLCLOCK]", self.wallClock ). \
                              replace("[MDRUNS]", self.mdRuns )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open("cases/"+caseName+"/submit_run.sh","w");        
        FILE.write( TEMP );
        FILE.close();

        # Show the file
        print TEMP
        
        # Let user review results
        self.confirmProgress()
        
    # Create all the amber input files for a case
    def amberCreateInput( self, caseName ):
        
        # User information
        self.printStage( "Stage 6, Case: "+caseName+". Creating Amber input files" )  

        # The template files for the amber imput files
        templateFiles = [
            self.PLMDHOME+"/src/templates/explicit_min.txt",
            self.PLMDHOME+"/src/templates/explicit_heat.txt",
            self.PLMDHOME+"/src/templates/explicit_equil.txt"
        ]
        
        
        # Set the QM region of this case
        self.calcQMregion( caseName )
        
        # Go through each template file
        for templateFile in templateFiles:

            # Load templates, change variables, and save in case folder
            TEMPLATE = open(templateFile, 'r')
            TEMP = TEMPLATE.read().replace("[NTC]", self.ntc ). \
                                  replace("[NTF]", self.ntf ). \
                                  replace("[QMCHARGE]", self.qmCharge ). \
                                  replace("[QMTHEORY]", self.qmTheory ). \
                                  replace("[QMREGION]", self.qmRegion ). \
                                  replace("[TIMESTEPS]", self.timestepNumber ). \
                                  replace("[DT]", str(self.timestepSize) ). \
                                  replace("[QMSHAKE]", self.qmShake ). \
                                  replace("[TIMESTEPPERFRAME]", str(self.timestepPerFrame) )
            TEMPLATE.close()
                                           
            # Save the input file with same name, but change extension to .in
            saveFile = os.path.basename(templateFile).split(".")[0]+".in"                                    
            FILE = open("cases/"+caseName+"/in_files/"+saveFile,"w");
            FILE.write( TEMP );
            FILE.close();
            
            # Show user the submission file
            print "\n"+TEMP+"\n"
            
            # Let user review results
            self.confirmProgress()
            

     # Function which analyses a final pdb file and figures out the QM region (ligand region)
    def calcQMregion( self, caseName ):
        
        # Open the pdb file created by LEaP
        with open("cases/"+caseName+"/pdb-files/finalLEaP_nowat.pdb",'r') as fl:
            pdb = fl.readlines()
        
        qmRegion = []
        
        # Go throug the file and find all residues having the resname of the ligand
        for line in pdb:
            if line[17:20] in self.ligandResnames:
                qmRegion.append( str(int(line[22:26])) )
        
        # Define the region string, as per Amber specifications
        if not qmRegion:
            
            # List was empty, not QM region
            self.qmRegion = ""
        
        else:
            
            # Set the QM region to the start-end ligand residues
            self.qmRegion = ":"+qmRegion[0]+"-"+qmRegion[ len(qmRegion)-1 ] 
            
        
         
    # Function to create a LEaP input file
    def leapCreateInput( self, caseName ):
        
        # User information
        self.printStage( "Stage 5, Case: "+caseName+". Creating input file for LEaP" )        
        
        # Forcefield loading
        ffString = "source "+self.AMBERHOME+"/dat/leap/cmd/"+self.ff 
        ffString += "\nsource "+self.AMBERHOME+"/dat/leap/cmd/leaprc.gaff"
        #ffString += "\nsource "+self.AMBERHOME+"/dat/leap/parm/frcmod.ionsjc_spce"

        # Structures Loading. Start with peptide
        structureString = "compound = loadpdb "+self.peptide
        
        # Go through the ligands
        for i in range(0,len(self.ligandsForLEaP)):
            
            # Import ligand
            structureString += "\nligand"+str(i)+" = loadmol2 "+self.ligand
            
            # Translate to correct position
            if self.noTranslate == False:
                structureString += "\ntranslate ligand"+str(i)+" {"+str(self.ligandsForLEaP[i][0])+" "+str(self.ligandsForLEaP[i][1])+" "+str(self.ligandsForLEaP[i][2])+" }"
            
            # Combine with previous structure
            structureString += "\ncompound = combine { compound ligand"+str(i)+" }"
        
        # Explicit solvent case: Equilibration
        TEMPLATE = open( self.PLMDHOME+"/src/templates/LEaP_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FORCEFIELDS]", ffString ). \
                              replace("[STRUCTURES_IMPORT]", structureString ). \
                              replace("[FOLDER]", "cases/"+caseName )
        TEMPLATE.close()
        FILE = open("cases/"+caseName+"/LEaP.ff","w");
        FILE.write( TEMP );
        FILE.close();        
        
        # Show user the submission file
        print TEMP
        
        # Let user review results
        self.confirmProgress()
        
         
    # Run leap input for case
    def leapRunInput( self, caseName ):
        
        # Run submission file, redirect to log file
        os.system(self.xLEaP + "cases/"+caseName+"/LEaP.ff &> cases/"+caseName+"/LEaP_setup.log" )
            
        # User information
        self.printStage( "Stage 5, Case: "+caseName+". LEaP Run Results" ) 
        with open("cases/"+caseName+"/LEaP_setup.log", 'r') as fin:
            print fin.read()
        
        # Let user review results
        self.confirmProgress()

    # Create main folder for the cases
    def createMainFolder( self ):
                
        # User information
        self.printStage( "Stage 2, Creating main cases/ folder" )
        
        # Check if folder already exists, confirm deletion
        self.createFolder( "cases" )        
            
    # Create simulation folder structure
    def creatCaseFolder(self, folderString):
        
        # User information
        self.printStage( "Stage 3, Case: "+folderString+". Setting up case folder structures" )
            
        # Setup required folders
        self.createFolder( "cases/"+folderString )            
        self.createFolder( "cases/"+folderString+"/md-files" )            
        self.createFolder( "cases/"+folderString+"/in_files" )            
        self.createFolder( "cases/"+folderString+"/md-logs" )            
        self.createFolder( "cases/"+folderString+"/pdb-files" )            
        
        # Let the user keep up
        if self.quiet == False:
            raw_input("Press any key to continue")            
            