#!/usr/bin/python
import os, re
import plmd 

class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Load the config file
        self.config = plmd.PLMD_Config( config ) 
        
        # Confirm with user
        self.printStage("Step 1: Starting up PLMD. Submission file parameters:")

        # Add GPU to nodecontrol if applicable
        if self.config.gpuEnabled == True:
            self.config.nodeControl += ":gpus="+str(self.config.gpuCores)
            
        # AMD setup variable
        self.aMDinput = ""
            

        print "\n== Submission Parameters"
        print "========================"
        print "submissionName: " + self.config.name
        print "nodeControl: " + self.config.nodeControl
        print "wallClock: " + self.config.wallClock
        print "mdRuns: " + self.config.mdRuns

        # Confirmation from user
        if self.config.quiet == False:
            var = raw_input("\nPlease confirm these submission parameters with any key press. Press 'n' to discontinue")
            if var == 'n':
                raise Exception('Submission was not confirmed')

    
    # Create submission file for submitting case to HPC queue
    def hpcCreateSubmission( self, caseName ):
        
        # Create in-files
        self.amberCreateInput( caseName )
        
        # User information
        self.printStage( "Stage 3, Case: "+caseName+". Creating HPC submission files" )          
        caseID = caseName.split("/")[-1] 
        
         # Create new submission file
        TEMPLATE = ""
        if self.config.gpuEnabled == True:
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/explicit_gpu_submit.txt", 'r')
        else:
            TEMPLATE = open( self.config.PLMDHOME+"/src/templates/explicit_submit.txt", 'r')
                        
        # What to call the logfile
        amdLogFile = ""
        if self.config.amdEnabled == True:
            amdLogFile = "outAMD"
        else:
            amdLogFile = "outMD"
        
        # Replace stuff within
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseName  ). \
                              replace("[NAME]", self.config.name+"_"+caseID  ). \
                              replace("[CPUCONTROL]", self.config.nodeControl ). \
                              replace("[WALLCLOCK]", self.config.wallClock ). \
                              replace("[MDRUNS]", self.config.mdRuns ). \
                              replace("[LOGFILENAME]", amdLogFile )
        TEMPLATE.close()
                              
        # Write the submission file
        FILE = open(caseName+"/submit_run.sh","w");        
        FILE.write( TEMP );
        FILE.close();
        
        print "Create submission file: "+caseName+"/submit_run.sh"
        
    
        
    # Create MMPBSA Submission file
    def hpcMMPBSASubmissionCreate(self, caseName, receptorCase, ligandCase ):
        
        # Create in-files
        self.amberCreateInput( caseName )
        
        # Run ante-MMPBSA.py
        self.runAnteMMPBSA(caseName)
        
        # User information
        self.printStage( "Stage 3, Case: "+caseName+". Creating HPC MMPBSA submission files. Receptor case: "+receptorCase )          
        caseID = caseName.split("/")[-1] 
                      
        # Add all trajectory files to ptraj script
        self.num_files = self.getNumberOfFiles( caseName+'/md-files/' ) 
        complexFiles = ""
        for i in range(1,self.num_files):
            complexFiles += caseName+'/md-files/equil'+ str(i)+ ".mdcrd "

        self.num_files = self.getNumberOfFiles( receptorCase+'/md-files/' ) 
        receptorFiles = ""
        for i in range(1,self.num_files):
            receptorFiles += receptorCase+'/md-files/equil'+ str(i)+ ".mdcrd "

        self.num_files = self.getNumberOfFiles( ligandCase+'/md-files/' ) 
        ligandFiles = ""
        for i in range(1,self.num_files):
            ligandFiles += receptorCase+'/md-files/equil'+ str(i)+ ".mdcrd "                       
            
        # Replace stuff within
        TEMPLATE = open( self.config.PLMDHOME+"/src/templates/mmpbsa_submit.txt", 'r')
        TEMP = TEMPLATE.read().replace("[FOLDER]", caseName  ). \
                              replace("[RECEPTORFOLDER]", receptorCase  ). \
                              replace("[LIGANDFOLDER]", ligandCase  ). \
                              replace("[NAME]", self.config.name+"_"+caseID  ). \
                              replace("[CPUCONTROL]", self.config.nodeControl ). \
                              replace("[WALLCLOCK]", self.config.wallClock ). \
                              replace("[CASEID]", str(caseName.split("/")[-1]) ). \
                              replace("[COMPLEXFILES]", complexFiles ). \
                              replace("[RECEPTORFILES]", receptorFiles ). \
                              replace("[LIGANDFILES]", ligandFiles )
        TEMPLATE.close()
                              
        # Create folder for this
        self.createFolder( caseName+"/mmpbsa" , True )                               
                              
        # Write the submission file
        FILE = open(caseName+"/submit_mmpbsa.sh","w");        
        FILE.write( TEMP );
        FILE.close();
        
        print "Create submission file: "+caseName+"/submit_mmpbsa.sh"

    # Get number of files
    def getNumberOfFiles( self, path ):        
        return len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and ".mdcrd" in f and "equil" in f] ) 
            
    
    def runAnteMMPBSA(self, caseName):

        # Delete old files?
        #os.system("rm -rf "+caseName+"/md-files/complex.prmtop "+caseName+"/md-files/receptor.prmtop "+caseName+"/md-files/ligand.prmtop")        
        
        command = "ante-MMPBSA.py \
        -p "+caseName+"/md-files/peptide.prmtop \
        -c "+caseName+"/md-files/complex.prmtop \
        -r "+caseName+"/md-files/receptor.prmtop \
        -l "+caseName+"/md-files/ligand.prmtop \
        -s \":WAT\" \
        -n \""+self.qmRegion+"\""
        os.system(command)
        
    # Submit to HPC cluster
    def hpcSubmission( self, caseName ):
        
        # User information
        self.printStage( "Stage 4, Case: "+caseName+". Submitting to HPC" )          
        
        # Do submission        
        os.system( "qsub "+caseName+"/submit_run.sh" )
        
    # Submit MMPBSA run
    def hpcMMPBSASubmission(self, caseName):
        
         # User information
        self.printStage( "Stage 4, MMPBSA Case: "+caseName+". Submitting to HPC" )          
        
        # Do submission        
        os.system( "qsub "+caseName+"/submit_mmpbsa.sh" )
             
    # Create all the amber input files for a case
    def amberCreateInput( self, caseName ):
        
        # User information
        self.printStage( "Stage 2, Case: "+caseName+". Creating Amber input files" )  

        # The template files for the amber imput files
        templateFiles = [
            self.config.PLMDHOME+"/src/templates/explicit_min.txt",
            self.config.PLMDHOME+"/src/templates/explicit_heat.txt",
            self.config.PLMDHOME+"/src/templates/explicit_equil.txt",
            self.config.PLMDHOME+"/src/templates/explicit_mmpbsa.txt"
        ]
        
        # Open the pdb file created by LEaP to find residues
        self.ligandResnames = []
        self.peptideResnames = []
        peptides = 0
        with open(caseName+"/pdb-files/finalLEaP_nowat.pdb",'r') as fl:
            for line in fl:
                
                # Check for TER commands
                if "TER" in line:
                    
                    # Increase count
                    peptides += 1
                else:
                    
                    # If we're not done with peptide, add resname                
                    if peptides < self.config.peptideCount:
                        if line[17:20] not in self.peptideResnames:
                            self.peptideResnames.append( line[17:20] )
                    else:
                        if line[17:20] not in self.ligandResnames:
                            self.ligandResnames.append( line[17:20] )
                            
                    
        
        # Set the QM region of this case
        self.calcQMregion( caseName )
        
        # Go through each template file
        for templateFile in templateFiles:

            # Enable quantum variable
            if self.config.ligandCount <= 0 or self.config.qmEnable == False:
                self.config.qmEnable = 0
            else:
                self.config.qmEnable = 1
                
            # Special things on equilfile
            if "equil" in templateFile:
                
                # GPU Optimization
                if self.config.gpuEnabled == True:
                    self.config.ntt = "3"
                    self.config.ntb = "1"
                    self.config.ntp = "0"
                    self.config.gamma_ln = "2.0"
                    
                # Enable aMD
                if self.config.amdEnabled == True:
                    
                    # Check for the latest equil log file to get aMD data. Otherwise abort
                    ePot, eDih = 0,0
                    for subdir,dirs,files in os.walk( caseName+"/md-logs/" ):
                        for filename in files:
                            if filename == "outMD1.log":
                                with open( subdir+filename , "r") as fi:
                                    startSearch = False
                                    for line in fi:
                                        
                                        # Start Search
                                        if "A V E R A G E S   O V E R" in line:
                                            startSearch = True
                                            
                                        # Do Search, only take first
                                        if startSearch == True:
                                             m1 = re.search('EPtot(\s+?)=(\s+?)(-?\d+\.?\d+)', line)
                                             if m1 and ePot == 0:
                                                ePot = float(m1.group(3))
                                             m2 = re.search('DIHED(\s+?)=(\s+?)(-?\d+\.?\d+)', line)
                                             if m2 and eDih == 0:
                                                eDih = float(m2.group(3))
                                            
                                        # End Search
                                        if "Density" in line:
                                            startSearch = False
                    
                    # Check that we found values
                    if ePot == 0 and eDih == 0:
                        raise Exception('To run aMD, a outMD1.log file must be present. This file is needed for information about the energies in the system.') 
                                
                    # Confirm aMD parameters
                    self.printStage( "Stage 2.5, Case: "+caseName+". aMD settings" )  
                    
                    # Get residues & atoms in the system
                    atoms, residues = 0,[]
                    with open( caseName+"/pdb-files/finalLEaP.pdb",'r' ) as fi:
                        for line in fi:
                            if "ATOM" in line:
                                atoms += 1
                            if line[17:20] in self.peptideResnames:
                                resID = str(int(line[22:26]))
                                if resID not in residues:
                                    residues.append( resID )
                        
                    # Input data
                    print "ATOMS: "+str(atoms)
                    print "RESIDUES: "+str(len(residues))
                    print "EPOT: "+str(ePot)
                    print "DIHED: "+str(eDih)
                    
                    # alphaP Calc
                    self.alphaP = self.config.ePA * atoms 
                    print "alphaP = "+str(self.config.ePA)+" * "+str(atoms)+" = "+str(self.alphaP)+" kcal mol-1"
                    
                    # EthreshP Calc
                    self.EthreshP = ePot + self.alphaP
                    print "EthreshP = "+str(ePot)+" + "+str(self.alphaP)+" = "+str(self.EthreshP)+" kcal mol-1"
                    
                    # EthreshD Calc
                    self.EthreshD = eDih + self.config.ePR * len(residues)
                    print "EthreshD = "+str(eDih)+" + "+str( self.config.ePR )+" * "+str(len(residues))+" = "+str(self.EthreshD)+" kcal mol-1"
                    
                    # alphaD Calc
                    self.alphaD = self.config.aDf * self.config.ePR * len(residues)
                    print "alphaD = "+str( self.config.aDf ) + " * " + str(self.config.ePR) + " * "+str(len(residues))+" = "+str(self.alphaD)+" kcal mol-1"

                    # Create entry for the input file
                    self.aMDinput = ",iamd="+str(self.config.iamd)+\
                               ",ethreshd="+str(self.EthreshD)+\
                               ",alphad="+str(self.alphaD)+\
                               ",ethreshp="+str(self.EthreshP)+\
                               ",alphap="+str(self.alphaP)
                    print self.aMDinput

                    # Confirm   
                    if self.config.quiet == False:                     
                        self.confirmProgress()            

            # Load templates, change variables, and save in case folder
            TEMPLATE = open(templateFile, 'r')
            TEMP = TEMPLATE.read().replace("[NTC]", self.config.ntc ). \
                                  replace("[NTF]", self.config.ntf ). \
                                  replace("[NTB]", self.config.ntb ). \
                                  replace("[NTT]", self.config.ntt ). \
                                  replace("[NTP]", self.config.ntp ). \
                                  replace("[GAMMALN]", self.config.gamma_ln ). \
                                  replace("[QMCHARGE]", self.config.qmCharge ). \
                                  replace("[QMTHEORY]", self.config.qmTheory ). \
                                  replace("[QMREGION]", self.qmRegion ). \
                                  replace("[TIMESTEPS]", self.config.timestepNumber ). \
                                  replace("[DT]", str(self.config.timestepSize) ). \
                                  replace("[PEPTIDERESI]", str(self.peptideRegion) ). \
                                  replace("[EABLEQM]", str(self.config.qmEnable) ). \
                                  replace("[QMSHAKE]", self.config.qmShake ). \
                                  replace("[AMDsetup]", self.aMDinput ). \
                                  replace("[COMPLEXIDS]", self.complexids ). \
                                  replace("[COMPLEXCHARGE]", str(self.config.qmCharge) ). \
                                  replace("[LIGANDCHARGE]", str(self.config.qmCharge) ). \
                                  replace("[INTERVAL]", str(self.config.mmpbsaInterval) ). \
                                  replace("[TIMESTEPPERFRAME]", str(self.config.timestepPerFrame) )
            
            
            # If not QM, delete qmmm dict from TEMP
            if self.config.qmEnable == 0:
                
                # Must be compiled first, so as to use DOTALL that will match newlines also
                TEMP = re.sub(re.compile('&qmmm(.+)\s/\n', re.DOTALL), "", TEMP )
                                           
            # Save the input file with same name, but change extension to .in
            saveFile = os.path.basename(templateFile).split(".")[0]+".in"                                    
            FILE = open(caseName+"/in_files/"+saveFile,"w");
            FILE.write( TEMP );
            FILE.close();
            
            
     # Function which analyses a final pdb file and figures out the QM region (ligand region)
    def calcQMregion( self, caseName ):
        
        # Open the pdb file created by LEaP
        with open(caseName+"/pdb-files/finalLEaP_nowat.pdb",'r') as fl:
            pdb = fl.readlines()
        
        # Set QM & peptide region
        qmRegion = []
        peptideRegion = []
        
        complexIDs = []
        receptorIDs = []
        ligandIDs = []
        
        # Go throug the file and find all residues having the resname of the ligand
        for line in pdb:
            if line[22:26]:
                resID = str(int(line[22:26]))
                if line[17:20] in self.ligandResnames:
                    qmRegion.append( resID )
                    if resID not in ligandIDs:
                        ligandIDs.append(resID)
                elif line[17:20] in self.peptideResnames:
                    peptideRegion.append( resID )
                    if resID not in receptorIDs:
                        receptorIDs.append(resID)
                if resID not in complexIDs:
                    complexIDs.append(resID)
        
        # Define the region string, as per Amber specifications
        if not qmRegion:
            
            # List was empty, not QM region
            self.qmRegion = ""
        
        else:
            
            # Set the QM region to the start-end ligand residues
            self.qmRegion = ":"+qmRegion[0]+"-"+qmRegion[ len(qmRegion)-1 ] 
            
        # Set the peptide region
        if peptideRegion:
            self.peptideRegion = ":"+peptideRegion[0]+"-"+peptideRegion[-1] 
        else:
            self.peptideRegion = ""
            
        # Get complex IDs for MMPBSA       
        self.complexids = ";".join(ligandIDs)
            