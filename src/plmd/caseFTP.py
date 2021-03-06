#!/usr/bin/python
import os
import tarfile
import plmd
import ftplib

# This is the overall Analysis class which merges trajectories,
# manages the analysis handler, and emails the final results
class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Save config parameters
        self.config = config 
    
    # Change directories
    def chdir(self, ftp, dir): 
        if self.ftp_hasDir(ftp, dir) is False: # (or negate, whatever you prefer for readability)
            ftp.mkd(dir)
        ftp.cwd(dir)    
    
    # Check if directory exists (in current location)
    def ftp_hasDir(self, ftp, dir):
        filelist = []
        ftp.retrlines('LIST',filelist.append)
        for f in filelist:
            if f.split()[-1] == dir and f.upper().startswith('D'):
                return True
        return False
        
    # Zip a given directory
    def zipDirectory( self, archive, folder ):
        self.printStage("Compressing "+folder + " to "+ archive + ".tar")
        tfile = tarfile.open(archive+".tar", "w:gz")
        tfile.add(folder)
        tfile.close()
               
    # Send a given file to the user ftp
    def shipTar( self, filepath, caseDir, postpend = "data" ):
        
        # User info
        serverFile = "tars/"+self.config.name+'_'+caseDir.split("/")[-1]+"_"+postpend+".tar"
        self.printStage("Sending file to FTP: "+filepath + " to " + serverFile)
        
        # Upload file
        session = ftplib.FTP( self.config.server , self.config.ftpUser , self.config.password )
        file = open(filepath,'r')                  
        session.storbinary( 'STOR ' + serverFile , file)    
        file.close()                                   
        session.quit()        
        
        # Done
        print "File Uploaded"

    # Send a directory to the user ftp
    def shipDir( self, shipDir , identifier = ""):
        
        # User info
        serverPath = self.config.name+'_'+identifier
        self.printStage("Sending directory to FTP: " + shipDir + " to " + serverPath)
        
        # Upload file
        session = ftplib.FTP( self.config.server , self.config.ftpUser , self.config.password )

        # Change the directory
        self.chdir( session, serverPath )
        
        # Go through the directories
        for subdir,dirs,files in os.walk( shipDir ):
            if files:
                for filename in files:
                    filepath = subdir + "/" + filename
                    print "Sending: ",filepath
                    file = open(filepath,'r')                  
                    session.storbinary( 'STOR ' + filename , file)    
                    file.close() 
                                  
        # Close the session
        session.quit()        
        