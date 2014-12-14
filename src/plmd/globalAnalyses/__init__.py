import MDAnalysis, os, re, sys
import plmd
import cluster, energy, pca, kld, endToEnd, RMSdFrequency, dihedral

# The analysis handler provides the interface to all the analysis modules
class analysisHandler (plmd.PLMD_module):
    
    def __init__(self, config, dataFiles, analysisDir ):

        # Save config parameters
        self.config = config
        
        # Directory for results
        self.resultDir = analysisDir

        # The directory we're working in 
        self.dataFiles = dataFiles

        # Load Trajectory of first case only. Just to get backbone structure
        self.mdTrajectories = []
        for caseDir in self.dataFiles['caseDirs']:
            self.mdTrajectories.append( MDAnalysis.Universe( 
                caseDir+"/md-files/peptide_nowat.prmtop", 
                caseDir+"/mergedResult.dcd") 
            ) 
        self.backbone = self.mdTrajectories[0].selectAtoms('protein and backbone')    
        
    # Run all the analyses modules
    def runAll( self ):
    
        # Do the PCA analysis
        try:
            pca.runAnalysis( self.dataFiles[ 'caseDirs' ] , self.resultDir, self.config.noReweight )
        except Exception as e:
            print "Failed pca analysis",e        
        
    
        # Cluster Comparisons
        try:
            self.printStage( "Comparing and plotting clusters")
            cluster.runAnalysis( self.dataFiles[ 'caseDirs' ] , self.resultDir, self.config.noReweight );
        except Exception as e:
            print "Failed cluster comparison & plotting",e    
            sys.exit(2)    
        
        # Plot Dihedral
        try:
            self.printStage( "Plotting dihedrals")
            dihedral.runAnalysis( self.dataFiles[ 'caseDirs' ] , self.resultDir, self.backbone, self.config.noReweight );
        except Exception as e:
            print "Failed dihedral analysis",e  
    
        # RMSd frequency
        try:
            RMSdFrequency.runAnalysis( self.dataFiles[ 'caseDirs' ] , self.resultDir, self.config.noReweight )    
        except Exception as e:
            print "Failed rmsd frequency analysis",e    
    
        # Run KLD analysis
        try:
            kld.runAnalysis( self.dataFiles['caseDirs'] , self.resultDir, self.mdTrajectories, self.backbone )    
        except Exception as e:
            print "Failed kld analysis",e
        
        
        # Plot Energies
        try:
            if self.config.noEnergy == False:
                self.printStage( "Plotting energies")
                energy.runAnalysis( "Kinetic Energies", self.dataFiles[ "summary.EKTOT" ] , self.resultDir );
                energy.runAnalysis( "Potential Energies", self.dataFiles[ "summary.EPTOT" ] , self.resultDir );
                energy.runAnalysis( "Total Energies", self.dataFiles[ "summary.ETOT" ] , self.resultDir );
        except Exception as e:
            print "Failed energy analysis",e
         
        # Plot all end-to-end distances on top of each other
        try:
            endToEnd.runAnalysis( self.dataFiles[ "dist_end_to_end.list.timeCorrected" ] , self.resultDir ) 
        except Exception as e:
            print "Failed endtoend analysis",e
        
    # Create and run ptraj file
    def loadPtrajTemplates( self ):

        # Get number of residues
        numOfResidues = self.backbone.numberOfResidues()

        # Load cpptraj template files for each case
        for caseDir in self.dataFiles[ 'caseDirs' ]:

            # Get all the global cpptraj templates
            globalTemplates = []
            for file in os.listdir(self.config.PLMDHOME+"/src/templates/"):
                m = re.search('cpptraj_analysis_global_(.+?).txt', file)
                if m:
                    globalTemplates.append( m.group(1) )

            # We run ptraj twice, once going through all frames (for analyses requiring that)
            # and once limited the trajectory to just enough points for a good plot
            for ptrajType in globalTemplates:        
                
                # Create new submission file
                TEMPLATE = open( self.config.PLMDHOME+"/src/templates/cpptraj_analysis_global_"+ptrajType+".txt", 'r')
                TEMP = TEMPLATE.read().replace("[FOLDER]", caseDir  ). \
                                       replace("[FIRSTRESI]", "1" ). \
                                       replace("[LASTRESI]", str(numOfResidues) ). \
                                       replace("[LASTID]", str( numOfResidues + self.config.ligandCount ) )
                TEMPLATE.close()
                                      
                # Write the submission file
                FILE = open(caseDir+"/ccptraj_analysis_"+ptrajType+".ptraj","w");        
                FILE.write( TEMP );
                FILE.close();
            
    def runPtrajCaseMerging( self ):
        
        # Combine all the trajectories
        trajectoryCombination = ""
        for caseDir in self.dataFiles['caseDirs']:
            trajectoryCombination += "\ntrajin "+caseDir+"/mergedResult.dcd 1 last 1"
        
        # Create cpptraj files. Do for all directories
        for caseDir in self.dataFiles['caseDirs']:
            TEMPLATE = open( caseDir+"/ccptraj_analysis_merge.ptraj", 'r')
            TEMP = TEMPLATE.read().replace("[INPUTTRAJ]", trajectoryCombination ). \
                                   replace("[OUTPUTDIR]", self.resultDir+"/data/") 
            TEMPLATE.close()
                                  
            # Write the submission file
            FILE = open(caseDir+"/ccptraj_analysis_merge.ptraj","w");        
            FILE.write( TEMP );
            FILE.close();
        
        # Run the cpptraj utility
        os.system( "$AMBERHOME/bin/cpptraj -p "+self.dataFiles['caseDirs'][0]+"/md-files/peptide_nowat.prmtop -i "+self.dataFiles['caseDirs'][0]+"/ccptraj_analysis_merge.ptraj" )
