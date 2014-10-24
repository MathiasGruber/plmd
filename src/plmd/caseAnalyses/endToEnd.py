#!/usr/bin/python
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor ):

     # User info
    print "Creating plot of end-to-end distance"    
    distanceFile = caseDir+"/analysis/data/dist_end_to_end.list"
    
    # Print time corrected files for later plotting

    with open( distanceFile , "r" ) as fl:
        lines = fl.readlines()
    with open( distanceFile+".timeCorrected" , "w") as fo:
        n = 0
        for line in lines:
            if n > 0:
                temp = line.split()
                fo.write( str(float(temp[0])*timeFactor) + "\t" + temp[1] +"\n")
            else:
                fo.write( line )
            n += 1
    

    # B factor plot
    myPlot.plotData( 
        caseDir+"/analysis/plots" , 
        "End to End Residue Distance", 
        ["End to End"], 
        [distanceFile] , 
        "Distance ($\AA$)", 
        xFactor = timeFactor,
        skipLines = 1 
    )
