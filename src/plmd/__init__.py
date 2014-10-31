import os,shutil
import plmd.defaultConfig

class PLMD_module:
        
    # Shows the user some text and prompts him to accept it to continue script
    def confirmProgress(self):
        
        # Let user review results
        var = raw_input("Please review above and confirm with any key press. Press 'n' to discontinue")
        if var == 'n':
            raise Exception('Progress was disconinued at the request of the user')
            
    # A function for printing the current stage of the process to the user
    def printStage( self,info ):
        print "\n"+"="*len(info)
        print info
        print "="*len(info)+"\n"
        
    # Function for creation of folders
    def createFolder( self, folderDir , silent = False ):
        
        # Check if the folder is already there
        if os.path.isdir( folderDir ) == True:
            if silent == False:
                var = raw_input("This will delete any previous data in the "+folderDir+"/ folder. Confirm with any key press. Press 'n' to discontinue")
                if var == 'n':
                    raise Exception('You opted not to delete the '+folderDir+'/ folder')
                
            # Delete all in old folder
            shutil.rmtree( folderDir )
        
        # Create new folder
        os.mkdir( folderDir )
        
        # Print info
        print "Created the folder: " + folderDir + "/"

# A class for holding the configuration settings
class PLMD_Config:
    
    # Initialize all configuration settings
    def __init__(self, configSettings):
        
        # First load default settings, then overwrite with parameter settings
        config = plmd.defaultConfig.getDefaultConfig()
        for section in configSettings.sections():
            for (key, value) in configSettings.items(section):
                config.set( section, key, configSettings.get( section, key ) )
                
        # PLMD settings
        self.name = config.get('plmd_settings', 'name')                
                
        # Save config parameters
        self.ligand = config.get('inputFiles', 'ligand')
        self.ligandCount = config.getint('inputFiles', 'ligandCount')
        self.pLigand = config.get('inputFiles', 'pLigand')
        self.peptide = config.get('inputFiles', 'peptide')
        self.pPeptide = config.get('inputFiles', 'pPeptide')
        self.peptideCount = config.getint('inputFiles', 'peptideCount')
        self.cases = config.getint('inputFiles', 'cases')
        self.noTranslate = config.getboolean('inputFiles', 'noTranslate')
        
        self.quiet = config.getboolean('inputFiles', 'quiet')
        
        self.ff = config.get('simulationParameters', 'forceField')
        self.qmEnable = config.getboolean('simulationParameters', 'qmEnable')
        self.qmCharge = config.get('simulationParameters', 'qmCharge')
        self.qmShake = config.get('simulationParameters', 'qmShake')
        self.qmTheory = config.get('simulationParameters', 'qmTheory')
        self.ntf = config.get('simulationParameters', 'ntf')
        self.ntc = config.get('simulationParameters', 'ntc')
        self.timestepPerFrame = config.getint('simulationParameters', 'timestepPerFrame')
        self.timestepSize = config.getfloat('simulationParameters', 'timestepSize')
        self.timestepNumber = config.get('simulationParameters', 'timestepNumber')
        
        self.nodeControl = config.get('submissionParameters', 'nodeControl')
        self.wallClock = config.get('submissionParameters', 'wallClock')
        self.mdRuns = config.get('submissionParameters', 'mdRuns')
        
        # GPU run
        self.gpuEnabled = config.getboolean('gpuRun', 'enable')       
        self.gpuCores =  config.getint('gpuRun', 'gpus')       
        
        # Analysis configuration
        self.noMerge = config.getboolean('analysisParameters', 'noMerge')
        self.noStrip = config.getboolean('analysisParameters', 'noStrip')
        self.noBlock = config.getboolean('analysisParameters', 'noBlock')
        self.noEnergy = config.getboolean('analysisParameters', 'noEnergy')
        self.noFTP = config.getboolean('analysisParameters', 'noFTP')
        self.clusterFrames = config.getfloat('analysisParameters', 'clusterFrames')
        
        # FTP configuration
        if config.has_option( "ftpConfiguration", "ftpPass" ):
            self.server = config.get('ftpConfiguration', 'server')
            self.ftpUser = config.get('ftpConfiguration', 'ftpUser')
            self.password = config.get('ftpConfiguration', 'ftpPass')
        else:
            self.noFTP = True
        
        # Get environment vars for PLMD & Amber. Will raise exceptions if not found
        self.PLMDHOME = os.environ["PLMDHOME"]
        self.AMBERHOME = os.environ["AMBERHOME"] 