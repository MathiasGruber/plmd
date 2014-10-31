#!/usr/bin/python
import os, shutil, re, sys
from copy import deepcopy
import structureManipulation, peptideConstructor, plmd 

class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Load the config file
        self.config = plmd.PLMD_Config( config )     
        
        # Variables for storing initial data about the components
        self.peptideCoordinates = []
        self.peptideCenterOfMass = []
        self.ligandCoordinates = []
        self.ligandCenterOfMass = []
        
        # Variables used during case setup
        self.ligandsForLEaP = []
        self.ligandResnames = []
        self.peptideResnames = []
        self.qmRegion = ""
        self.peptideRegion = ""
        
        # Variable for storing possible leap modifications
        self.leapMods = False
        
        # Ensure integrity of input files
        self.setupInputFiles()
        
        # Setup quick links to tools in amber
        self.xLEaP = self.config.AMBERHOME+"/bin/tleap -f "

        # Confirm with user
        self.printStage("Step 1: Starting up PLMD. Configuration file:")

        print "\n== Input Files"
        print "=============="
        print "ligand: " + self.config.ligand
        print "ligand count: " + str(self.config.ligandCount)
        print "peptide: " + self.config.peptide
        print "peptide count: " + str(self.config.peptideCount)
        print "cases: " + str(self.config.cases)

        print "\n== Simulation Parameters"
        print "========================"
        print "Forcefield: " + self.config.ff
        print "qmCharge: " + self.config.qmCharge 
        print "qmShake: " + self.config.qmShake
        print "qmTheory: " + self.config.qmTheory
        print "ntf: " + self.config.ntf 
        print "ntc: " + self.config.ntc
        print "timestepSize: " + str(self.config.timestepSize)
        print "timestepNumber: " + self.config.timestepNumber

        print "\n== Submission Parameters"
        print "========================"
        print "nodeControl: " + self.config.nodeControl
        print "wallClock: " + self.config.wallClock
        print "mdRuns: " + self.config.mdRuns
        
        if self.config.quiet == False:
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
        if self.config.pPeptide != None and self.config.pPeptide != "false":

            # User info
            print "Now creating peptide with LEaP"
            
            # Create & save the peptide
            peptideCreator = peptideConstructor.Creator()
            peptideCreator.createPeptide( self.config.pPeptide, self.config.peptideCount, "predefinedInput/peptides.pdb" )
            
            # Overwrite default configs
            self.config.peptide = "predefinedInput/peptides.pdb"
              
        # Check if we should use predefined ion from package
        if self.config.pLigand != None and self.config.pLigand != "false":
            
            # User info
            print "Now importing the predefined ion to be user to library."
            
            # Check that the file exists in pre-defined library
            if os.path.isfile( self.config.PLMDHOME+"/src/ligands/"+self.config.pLigand+"/"+self.config.pLigand+".mol2" ):
                
                # Copy and set
                shutil.copy( self.config.PLMDHOME+"/src/ligands/"+self.config.pLigand+"/"+self.config.pLigand+".mol2", "predefinedInput/"+self.config.pLigand+".mol2" )
                self.config.ligand = "predefinedInput/"+self.config.pLigand+".mol2"
                
                # Check if there are LEaP modifications to be loaded
                if os.path.isfile( self.config.PLMDHOME+"/src/ligands/"+self.config.pLigand+"/leapLigand.ff" ):
                    
                    # User info
                    print "Found a LEaP modification for this file. Saving for future reference"                    
                    
                    # Save the leap modifications
                    with open( self.config.PLMDHOME+"/src/ligands/"+self.config.pLigand+"/leapLigand.ff" ) as fl:
                        self.leapMods = fl.readlines()
                
            else:
                
                # Raise error
                raise Exception("Tried to include invalid pre-defined ligand.")                

    # Main functon for setting up all the cases
    def setupCases(self):
        
        # Get the coordinates, resnames and center of masses of peptide. Passed by reference.
        if self.config.peptideCount > 0:
            structureManipulation.processStructureFile( self.config.peptide , "pdb", self.peptideCoordinates, self.peptideResnames )
            self.peptideCenterOfMass = structureManipulation.centerOfMass( self.peptideCoordinates )
        else:
            raise Exception("A simulation without protein is senseless. If you are trying to debug the ion, please run Amber manually for increased control.") 
        
        # Get the coordinates, resnames and center of masses of ligand. Passed by reference.
        if self.config.ligandCount > 0:
            structureManipulation.processStructureFile( self.config.ligand , "mol2", self.ligandCoordinates, self.ligandResnames )
            self.ligandCenterOfMass = structureManipulation.centerOfMass( self.ligandCoordinates )
        
        
        # Create main folder for project
        self.createMainFolder()        
        
        # Go through all the cases to setup
        for i in range(0,self.config.cases):
            
            # Create folder structure
            self.creatCaseFolder( str(i) )            
            
            # User information
            self.printStage( "Stage 4, Case: "+str(i)+". Calculation ion translation vector" )
            
            # Calculate a random translation for the ligand (pass by value)
            # LEaP performs the translation using this vector
            if self.config.ligandCount > 0:
                self.ligandsForLEaP[:] = []
                for x in range(0, self.config.ligandCount):
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
        
    # Create all the amber input files for a case
    def amberCreateInput( self, caseName ):
        
        # User information
        self.printStage( "Stage 6, Case: "+caseName+". Creating Amber input files" )  

        # The template files for the amber imput files
        templateFiles = [
            self.config.PLMDHOME+"/src/templates/explicit_min.txt",
            self.config.PLMDHOME+"/src/templates/explicit_heat.txt",
            self.config.PLMDHOME+"/src/templates/explicit_equil.txt",
            self.config.PLMDHOME+"/src/templates/explicit_gpu_equil.txt"
        ]
        
        
        # Set the QM region of this case
        self.calcQMregion( caseName )
        
        # Go through each template file
        for templateFile in templateFiles:

            # Enable quantum variable
            if self.config.ligandCount <= 0 or self.config.qmEnable == False:
                self.config.qmEnable = 0
            else:
                self.config.qmEnable = 1

            # Load templates, change variables, and save in case folder
            TEMPLATE = open(templateFile, 'r')
            TEMP = TEMPLATE.read().replace("[NTC]", self.config.ntc ). \
                                  replace("[NTF]", self.config.ntf ). \
                                  replace("[QMCHARGE]", self.config.qmCharge ). \
                                  replace("[QMTHEORY]", self.config.qmTheory ). \
                                  replace("[QMREGION]", self.qmRegion ). \
                                  replace("[TIMESTEPS]", self.config.timestepNumber ). \
                                  replace("[DT]", str(self.config.timestepSize) ). \
                                  replace("[PEPTIDERESI]", str(self.peptideRegion) ). \
                                  replace("[EABLEQM]", str(self.config.qmEnable) ). \
                                  replace("[QMSHAKE]", self.config.qmShake ). \
                                  replace("[TIMESTEPPERFRAME]", str(self.config.timestepPerFrame) )
            TEMPLATE.close()
            
            # If not QM, delete qmmm dict from TEMP
            if self.config.qmEnable == 0:
                
                # Must be compiled first, so as to use DOTALL that will match newlines also
                TEMP = re.sub(re.compile('&qmmm(.+)\s/\n', re.DOTALL), "", TEMP )
                                           
            # Save the input file with same name, but change extension to .in
            saveFile = os.path.basename(templateFile).split(".")[0]+".in"                                    
            FILE = open("cases/"+caseName+"/in_files/"+saveFile,"w");
            FILE.write( TEMP );
            FILE.close();
            
            # Show user the submission file
            print "\n"+TEMP+"\n"
            
            # Let user review results
            if self.config.quiet == False:
                self.confirmProgress()
            

     # Function which analyses a final pdb file and figures out the QM region (ligand region)
    def calcQMregion( self, caseName ):
        
        # Open the pdb file created by LEaP
        with open("cases/"+caseName+"/pdb-files/finalLEaP_nowat.pdb",'r') as fl:
            pdb = fl.readlines()
        
        qmRegion = []
        peptideRegion = []
        
        # Go throug the file and find all residues having the resname of the ligand
        for line in pdb:
            if line[17:20] in self.ligandResnames:
                qmRegion.append( str(int(line[22:26])) )
            elif line[17:20] in self.peptideResnames:
                peptideRegion.append( str(int(line[22:26])) )
        
        # Define the region string, as per Amber specifications
        if not qmRegion:
            
            # List was empty, not QM region
            self.qmRegion = ""
        
        else:
            
            # Set the QM region to the start-end ligand residues
            self.qmRegion = ":"+qmRegion[0]+"-"+qmRegion[ len(qmRegion)-1 ] 
            
        # Set the peptide region
        self.peptideRegion = ":"+peptideRegion[0]+"-"+peptideRegion[-1] 
            
        
         
    # Function to create a LEaP input file
    def leapCreateInput( self, caseName ):
        
        # User information
        self.printStage( "Stage 5, Case: "+caseName+". Creating input file for LEaP" )        
        
        # Forcefield loading
        ffString = "source "+self.config.AMBERHOME+"/dat/leap/cmd/"+self.config.ff 
        ffString += "\nsource "+self.config.AMBERHOME+"/dat/leap/cmd/leaprc.gaff"
        
        # Any additional 
        if self.leapMods != False:
            for line in self.leapMods:
                ffString += line

        # Structures Loading. Start with peptide
        structureString = "compound = loadpdb "+self.config.peptide
        
        # Go through the ligands
        if self.config.ligandCount > 0:
            for i in range(0,len(self.ligandsForLEaP)):
                
                # Import ligand
                structureString += "\nligand"+str(i)+" = loadmol2 "+self.config.ligand
                
                # Translate to correct position
                if self.config.noTranslate == False:
                    structureString += "\ntranslate ligand"+str(i)+" {"+str(self.ligandsForLEaP[i][0])+" "+str(self.ligandsForLEaP[i][1])+" "+str(self.ligandsForLEaP[i][2])+" }"
                
                # Combine with previous structure
                structureString += "\ncompound = combine { compound ligand"+str(i)+" }"
        
        # Explicit solvent case: Equilibration
        TEMPLATE = open( self.config.PLMDHOME+"/src/templates/LEaP_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FORCEFIELDS]", ffString ). \
                              replace("[PLMDHOME]", self.config.PLMDHOME ). \
                              replace("[STRUCTURES_IMPORT]", structureString ). \
                              replace("[FOLDER]", "cases/"+caseName )
        TEMPLATE.close()
        FILE = open("cases/"+caseName+"/LEaP.ff","w");
        FILE.write( TEMP );
        FILE.close();        
        
        # Show user the submission file
        print TEMP
        
        # Let user review results
        if self.config.quiet == False:
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
        if self.config.quiet == False:
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
        if self.config.quiet == False:
            raw_input("Press any key to continue")            
            