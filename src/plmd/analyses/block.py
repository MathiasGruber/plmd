#!/usr/bin/python
import numpy, matplotlib
from pylab import plot,errorbar, subplot, xlabel, ylabel, savefig, plt
import MDAnalysis

# == BLOCK AVERAGING MODULE
# == For details on the mathematics behind block averaging, read:
# == Quantifying uncertainty and sampling quality in biomolecular simulations
# == Annu Rep Comput Chem. 2009 January 1; 5: 23-48
# =============================================

# Block Analysis
def blocked(universe, nblocks, analyze):
    size = universe.trajectory.numframes/nblocks
    blocks = []
    for block in xrange(nblocks):
        a = []
        for ts in universe.trajectory[block*size:(block+1)*size]:
            a.append(analyze(universe))
        blocks.append(numpy.average(a))
    blockaverage = numpy.average(blocks)
    blockstd = numpy.std(blocks)
    
    return nblocks, size, blockaverage, blockstd

# Function to use in block analysis
def rgyr(universe):
    return universe.selectAtoms('protein').radiusOfGyration()

# Function for running the actual analysis
def runAnalysis( saveDir, mdTrajectory , timePerFrame ):

    # Pure radius of gyration
    print "Calculating radius of gyration for trajectory"
    tmp = []
    for ts in mdTrajectory.trajectory:
        time = ts.frame*timePerFrame
        tmp.append([time,rgyr(mdTrajectory)])
    rg = numpy.array(tmp)

    # Run the block analysis    
    results = []
    for nblocks in xrange(5,20):
        print "Running for "+str(nblocks)+" blocks"
        results.append(blocked(mdTrajectory, nblocks, rgyr))
        r = numpy.array(results)

    # Print block analysis
    have_matplotlib = True
    if have_matplotlib:
        
        font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}
        plt.rc('font', **font)
        
        subplot(221)
        errorbar(r[:,0], r[:,2], yerr=r[:,3])
        xlabel("number of blocks")
        ylabel(r"$\langle R_{\rm{gyr}} \rangle$ ($\AA$)")
        
        subplot(222)
        errorbar(r[:,1], r[:,2], yerr=r[:,3])
        xlabel("block size")
        ylabel(r"$\langle R_{\rm{gyr}} \rangle$ ($\AA$)")

        subplot(223)
        plot(rg[:,0], rg[:,1])
        xlabel("Time")
        ylabel(r"$R_{\rm{gyr}}$ ($\AA$)")
        
        subplot(224)
        plot(r[:,1]*timePerFrame*1000, numpy.divide(r[:,3] , numpy.sqrt(r[:,0]) ) )
        xlabel("block length (ps)")
        ylabel("Blocked Standard Error")
        
        savefig(saveDir+"/analysis/plots/blocks.pdf")
        
        print "Wrote ./figures/blocks.{pdf,png}" % vars()