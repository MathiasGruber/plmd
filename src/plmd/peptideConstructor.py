import os

# Class that handles creation of new linear peptides using LEaP
class Creator:
    
    # Constructor for sanity checks
    def __init__( self ):
        
        # Make sure the required libraries are isntalled
        self.AMBERHOME = os.environ["AMBERHOME"] 
        self.PLMDHOME = os.environ["PLMDHOME"]
        self.xLEaP = self.AMBERHOME+"/bin/tleap -f "
        
    # Creates a linear peptide from sequence using LEaP
    def createPeptide( self , sequenceString, peptideCount, outFile ):
        
        # Check that output path exists
        outDir = os.path.dirname( outFile )
        if os.path.isdir( outDir ) or not outDir:
            
            # If we are in current directory
            if not outDir:
                outDir = "./"
            
            # Available residues       
            resNames = [ "ALA","ARG","ASN","ASX","CYS","GLU","GLN","GLX","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL" ]
            cTerm = ["C"+b for b in resNames]
            nTerm = ["N"+b for b in resNames]
    
            # Validate sequence            
            sequence = sequenceString.split()
            for entry in sequence:
                if entry not in resNames and entry not in cTerm and entry not in nTerm:
                    raise Exception("Invalid sequance supplied with residue entry: "+entry)
                    
            # Create input file
            ffString = "source "+self.AMBERHOME+"/dat/leap/cmd/leaprc.ff12SB"
    
            # Structures Loading. Start with initial single peptide
            structureString = "compound = sequence { "+sequenceString+" }"    
            if peptideCount > 1:
                for i in range( 0, peptideCount ):
                    structureString += "\npeptide"+str(i)+" = sequence { "+sequenceString+" }"    
                    structureString += "\ntranslate peptide"+str(i)+" {0 "+str(i*10)+" 0}"
                    structureString += "\ncompound = combine { compound peptide"+str(i)+" }"
                
            # Explicit solvent case: Equilibration
            TEMPLATE = open( self.PLMDHOME+"/src/templates/LEaP_peptide.txt", 'r')
            TEMP = TEMPLATE.read().replace("[FORCEFIELDS]", ffString ). \
                                  replace("[STRUCTURES_IMPORT]", structureString ). \
                                  replace("[OUTPUT]", outFile )
            TEMPLATE.close()
            FILE = open(outDir+"/LEaP.ff","w");
            FILE.write( TEMP );
            FILE.close();        
            
            # Show user the submission file
            # Run submission file, redirect to log file
            os.system(self.xLEaP + outDir+"/LEaP.ff &> "+outDir+"/LEaP_peptideCreation.log" )
                
            # User information
            print "Script finished. \n Check the file: "+outFile+" to inspect final result.  \n Check the file: "+outDir+"/LEaP_peptideCreation.log for debug information"
            
        else:
            
            raise Exception("Specified output directory does not exist")