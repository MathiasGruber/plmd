#!/usr/bin/python
import os
import numpy as np
import plmd.plotData as myPlot

# == Takes energy data from MD logs and plots it
# ==============================================

# Function for running the actual analysis
def runAnalysis( caseDir , timeFactor, noReweight = False ):

    # Check if the file exists
    if os.path.exists( caseDir+"/mmpbsa/energyPerFrame.dat" ):
        
        # User info
        print "Found MMPBSA Files. Plotting binding energy vs. time"
        
        # Recreate file to frame & energy
        energies = []
        startCounting = 0
        lines = 0
        with open(caseDir+"/mmpbsa/energyPerFrame.dat", "r") as fi:
            with open(caseDir+"/mmpbsa/energyPerFrame_cleaned.dat", "w") as fo:
                next(fi)
                next(fi)
                next(fi)
                for line in fi:
                    temp = line.split(",")
                    if line.strip() and temp[-1] and startCounting == 1:
                        energies.append(float(temp[-1]))
                        lines += 1
                        fo.write( str(lines) + "\t" + temp[-1] + "\n" )
                    if "DELTA Energy Terms" in line:
                        startCounting = 1
                        next(fi)

        # Common Plot command
        myPlot.plotData( 
            caseDir+"/analysis/plots" , 
            "MMPBSA Energy", 
            ["MMPBSA Energy"],
            [ caseDir+"/mmpbsa/energyPerFrame_cleaned.dat" ] , 
            "Energy [kcal / mol]", 
            xFactor = timeFactor,
            scatter = True ,
            legendFrame = 1,
            legendAlpha = 1
        )
 
       
        