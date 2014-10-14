#!/usr/bin/python
import os
import tarfile
import plmd

# Sending email modules
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders

# This is the overall Analysis class which merges trajectories,
# manages the analysis handler, and emails the final results
class Setup (plmd.PLMD_module):

    def __init__(self, config):
        
        # Save config parameters
        self.load_config( config )  
        
    # Zip a given directory
    def zipDirectory( self, archive, folder ):
        self.printStage("Compressing: "+folder)
        tfile = tarfile.open(archive+".tar", "w:gz")
        tfile.add(folder)
        tfile.close()
               
    # Email a given file to the user
    def emailFile( self, filepath, caseDir ):
        self.printStage("Emailing the file: "+filepath)
        
        # Message to be sent
        msg = MIMEMultipart()
        msg['Subject'] = 'PLMD: '+self.name+" : "+caseDir
        msg['From'] = self.fromEmail
        msg['To'] = self.toEmail  
        msg.attach( MIMEText("""
Results are in for: 
Simulation Name: """+self.name+"""
Case Folder: """+caseDir+"""
""") )        
        
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(filepath,"rb").read() )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(filepath))
        msg.attach(part)
        
        # Info to user
        print "Sending email to: "+self.toEmail
        print "Using smtp server:"+self.smtp
        
        # Connect to SMTP
        server = smtplib.SMTP( self.smtp , self.port )
        server.ehlo()
        server.starttls()
        server.login( self.fromEmail, self.password)
        
        # Ship off email
        server.sendmail(self.fromEmail, [self.toEmail], msg.as_string() )
        server.quit()
        
        # Done
        print "Email sent"
