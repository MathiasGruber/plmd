#!/usr/bin/python
import numpy, matplotlib
from pylab import plot,errorbar, subplot, xlabel, ylabel, savefig, plt
import MDAnalysis

# Function for running the actual analysis
def runAnalysis( saveDir, mdTrajectory , timePerFrame ):

    print "Doing some H-bond stuff"