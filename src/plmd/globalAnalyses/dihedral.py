#!/usr/bin/python
import plmd.plotData as myPlot
import re

# Do the plotting   
def doDihedralPlot( plotTitle, datafiles ):
    
    print "Creating plot of: "+plotTitle    
    
    myPlot.plotData( 
        "globalAnalyses" , 
        plotTitle , 
        datafiles['caseLabels'], 
        datafiles['filepaths'] , 
        "E (kcal/mol)" 
    )

# Function for running the actual analysis
def runAnalysis( dataFiles, backbone ):
    
    resNames = backbone.resnames()
    for k,v in dataFiles.items():
        if ("psi" in k or "phi" in k) and "timeCorrected" in k:
            resi = int(re.findall(r'\d+', k )[0])
            eTitle = ". Res ID: "+str(resi)+", Residue name: "+str(resNames[resi])
            doDihedralPlot( k.split()[-1]+" Angles" + eTitle, dataFiles[ k ] )   
    
    
 