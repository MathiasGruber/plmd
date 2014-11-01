#!/usr/bin/python
import os, sys
import plmd.generalAnalyses.componentAnalysis as pcaFuncs
import plmd.plotData as myPlot

# Function for running the actual analysis
def runAnalysis( caseDirs, resultsDir, trajectories, backbone ):
    
    # User info
    print "Doing KLD analysis."
    
    # RUN KLD & PCA ON TOTAL    
    
    # Use the cpptraj script in the first case
    cppTrajScript = caseDirs[0]+"/ccptraj_analysis_kld.ptraj"
    
    # Create KLD projections
    crdActions = ""
    caseNames = []
    i = 0
    startFrame = 1
    numOfResidues = backbone.numberOfResidues()
    for caseDir in caseDirs:
        frames = trajectories[i].trajectory.numframes
        lastFrame = startFrame + frames
        caseName = caseDir.split("/")[-1]
        caseNames.append(caseName)
        crdActions += "\ncrdaction crd1 projection T"+caseName+" modes "+resultsDir+"/data/evecs-ca.dat beg 1 end 20 :1-"+str(numOfResidues)+"&!@H= crdframes "+str(startFrame)+","+str(lastFrame)+" out "+resultsDir+"/data/T"+caseName+".dat" 
        startFrame = lastFrame      
        i += 1
        
    # Calculate Kullback-Leibler Divergence vs time for PC histograms 
    # Trajectories i and j, modes 1-5
    crdActions += "\n"
    
    # How many modes to analyze with KLD
    modesToAnalyze = 3   
    kldData = {}
    for i in range( 1,modesToAnalyze+1 ):
        kldData[ i ] = {"files":[],"labels":[]}
    
    # Create KLD part of cpptraj
    trash = []
    for firstCase in caseNames:
        trash.append(firstCase)
        for secondCase in caseNames:
            if secondCase not in trash:
                for i in range(1,modesToAnalyze+1):
                    
                    # Cpptraj Entry
                    outFile = resultsDir+"/data/KL-PC-M"+str(i)+"-"+firstCase+"-"+secondCase+".agr"
                    crdActions += "\nkde T"+firstCase+":"+str(i)+" kldiv T"+secondCase+":"+str(i)+" klout "+outFile+" bins 400 name KLD-M"+str(i)+"-"+firstCase+"-"+secondCase
                
                    # Save for plotting
                    kldData[i]['files'].append( outFile )   
                    kldData[i]['labels'].append( "Cases "+firstCase+" and "+secondCase )   
        
    # Calculate PC histograms
    histogramLabels, kdeHist, pcaHist = [],[],[]
    for trajname in caseNames:
        
        # Labels
        histogramLabels.append( "Case: "+trajname )        
        
        # Calculate PC histogram with KDE, trajectory 1, mode 1    
        crdActions += "\nkde T"+trajname+":1 out "+resultsDir+"/data/kde-PC-T"+trajname+".agr bins 200 name KDE1-T"+trajname
        kdeHist.append( resultsDir+"/data/kde-PC-T"+trajname+".agr" )
    
        # Calculate PC histogram, trajectory 1, mode 1
        crdActions += "\nhist T"+trajname+":1,*,*,*,200 out "+resultsDir+"/data/pca.hist.T"+trajname+".agr normint name HIST-1-T"+trajname 
        pcaHist.append( resultsDir+"/data/pca.hist.T"+trajname+".agr" )
    
    # Open the cpptraj script and add information
    TEMPLATE = open( cppTrajScript, 'r')
    TEMP = TEMPLATE.read().replace("[ANALYSISDIR]", resultsDir ). \
                           replace("[KLDINFO]", crdActions )
    TEMPLATE.close()
                          
    # Write the submission file
    FILE = open(cppTrajScript,"w");        
    FILE.write( TEMP );
    FILE.close();

    # Run the cpptraj utility
    os.system( "$AMBERHOME/bin/cpptraj -p "+caseDirs[0]+"/md-files/peptide_nowat.prmtop -i "+cppTrajScript )

    # Create PC histogram plots
    myPlot.plotData( 
        resultsDir+"/plots" , 
        "PC histograms from KDE", 
        histogramLabels , 
        kdeHist , 
        "Frequency (%)",
        scatter = True,
        skipLines = 8,
        xUnit = "PC1"
    )
    
    # Create PC histogram plots
    myPlot.plotData( 
        resultsDir+"/plots" , 
        "PC histograms from PCA", 
        histogramLabels , 
        pcaHist , 
        "Frequency (%)",
        scatter = True,
        skipLines = 8,
        xUnit = "PC1"
    )

    # Plot KLDs for each mode
    for i in range( 1,modesToAnalyze+1 ):
        myPlot.plotData( 
            resultsDir+"/plots" , 
            "KLD for mode "+str(i), 
            kldData[i]['labels'] , 
            kldData[i]['files'] , 
            "KLD",
            scatter = True,
            skipLines = 8,
            tightYlimits = True
        )


    # CREATE PCA PLOTS

    # Create eigenvalue file
    with open( resultsDir+"/data/evecs-cleaned.dat", "w" ) as fo:
        with open( resultsDir+"/data/evecs-ca.dat","r" ) as f:
            lookForVec = 1
            for line in f:
                temp = line.split()
                if len(temp) == 2 and temp[0] == str(lookForVec):
                    if lookForVec < 15:
                        fo.write( str(lookForVec)+"\t"+temp[1]+"\n" )
                        lookForVec += 1

    # Do the plotting
    pcaHandler = pcaFuncs.PCA( 
        resultsDir+"/plots/Combined_PCA.pdf"
    )

    # Do the plots of energy landscapes & distributions
    pcaHandler.plotPCA( 
        "Combined PCA Plot",              # Plot Title
        resultsDir+"/data/" ,             # Data Dir
        "pca12-ca",                       # Eigenvector file
        "evecs-cleaned.dat"               # Eigenvalue file
    )
    
    # Save the plot
    pcaHandler.savePlot()
    
    # Save an eigenvalue vs eigenvector plot
    print "Creating plot of eigenvalues"
    myPlot.plotData( 
        resultsDir+"/plots" , 
        "Eigenvalue Plot", 
        ["Eigenvalue"], 
        [resultsDir+"/data/evecs-cleaned.dat"] , 
        "AU",
        scatter = True
    )
    
    