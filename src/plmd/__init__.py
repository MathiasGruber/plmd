import os,shutil

class PLMD_module:
    
    def load_config(self, config):
    
        print"__INIT__ FILE LOADED"
        # Save config parameters
        self.ligand = config.get('inputFiles', 'ligand')
        self.ligandCount = config.getint('inputFiles', 'ligandCount')
        self.pLigand = config.get('inputFiles', 'pLigand')
        self.peptide = config.get('inputFiles', 'peptide')
        self.pPeptide = config.get('inputFiles', 'pPeptide')
        self.peptideCount = config.getint('inputFiles', 'peptideCount')
        self.cases = config.getint('inputFiles', 'cases')
        
        self.quiet = config.getint('inputFiles', 'quiet')
        
        self.ff = config.get('simulationParameters', 'forceField')
        self.qmCharge = config.get('simulationParameters', 'qmCharge')
        self.qmShake = config.get('simulationParameters', 'qmShake')
        self.qmTheory = config.get('simulationParameters', 'qmTheory')
        self.ntf = config.get('simulationParameters', 'ntf')
        self.ntc = config.get('simulationParameters', 'ntc')
        self.timestepSize = config.get('simulationParameters', 'timestepSize')
        self.timestepNumber = config.get('simulationParameters', 'timestepNumber')
        
        self.nodeControl = config.get('submissionParameters', 'nodeControl')
        self.wallClock = config.get('submissionParameters', 'wallClock')
        self.mdRuns = config.get('submissionParameters', 'mdRuns')
        
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
        
        # Get environment vars for PLMD & Amber. Will raise exceptions if not found
        self.PLMDHOME = os.environ["PLMDHOME"]
        self.AMBERHOME = os.environ["AMBERHOME"] 
        
        # Shows the user some text and prompts him to accept it to continue script
    def confirmProgress(self):
        
        # Let user review results
        if self.quiet == False:
            var = raw_input("Please review above and confirm with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Progress was disconinued at the request of the user')
            
    # A function for printing the current stage of the process to the user
    def printStage( self,info ):
        print "\n"+"="*len(info)
        print info
        print "="*len(info)+"\n"
        
    # Function for creation of folders
    def createFolder( self, folderDir ):
        
        # Check if the folder is already there
        if os.path.isdir( folderDir ) == True:
            var = raw_input("This will delete any previous data in the "+folderDir+"/ folder. Confirm with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('You opted not to delete the '+folderDir+'/ folder')
                
            # Delete all in old folder
            shutil.rmtree( folderDir )
        
        # Create new folder
        os.mkdir( folderDir )
        
        # Print info
        print "Created the folder: " + folderDir + "/"
