#!/usr/bin/env python
import math
import numpy as np
from random import gauss

# Function to make a random vector in n-dimensional space
def make_rand_vector(dims,length):
    vec = [gauss(0, 1) for i in range(dims)]
    mag = sum(x**2 for x in vec) ** .5 * (1/float(length))
    return [x/mag for x in vec]

# Read the coordinate list of a pdb or mol2 file
# Returns coordinates and also resnames of the file
# File specifications:
# PDB: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
# Mol2: http://www.tripos.com/data/support/mol2.pdf
def processStructureFile( filename , filetype, coor, resnames ):
    fl = open(filename,'r')
    pdb = fl.readlines()
    
    # Handle different fileformats
    if filetype == "pdb":
        i = 1
        for line in pdb:
            if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
                coor[i] = [ float(line[31:38]),float(line[39:46]),float(line[47:54]) ]
                resnames.append( line[17:20] )
                i = i + 1
    elif filetype == "mol2":
        i = 1
        readCoords = False
        for line in pdb:
            if "@" in line:
                readCoords = False
            if '@<TRIPOS>ATOM' in line :
                readCoords = True
                continue
            if readCoords == True:
                data = line.split()
                coor[i] = [ float(data[2]),float(data[3]),float(data[4]) ]
                resnames.append( data[7] )
                i = i + 1
    else:
        raise Exception('File '+filename+" could not be read") 
    return 1
    
# Calculate the center of mass (CoM) of a coordinate list (not accounting for density diffs)
def centerOfMass( structureCoords ):
    avgX,avgY,avgZ = np.array([]),np.array([]),np.array([])
    for k, strucCor in structureCoords.items():
        avgX = np.append( avgX, strucCor[0] )
        avgY = np.append( avgY, strucCor[1] )
        avgZ = np.append( avgZ, strucCor[2] )
    return [ np.average(avgX) , np.average(avgY), np.average(avgZ)]

# Calculate the distance between two points
def distance( point1, point2 ):
    return math.sqrt( (point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2 )

# Calculate vector between points
def calcVector( structureFrom, structureTo ):
    return [ structureFrom[0]-structureTo[0] , structureFrom[1]-structureTo[1] , structureFrom[2]-structureTo[2] ]

# Translate a structure (pass by reference)
def translateStructure( structureCoords, translateVector ):
    for i in structureCoords:
        structureCoords[i][0] = structureCoords[i][0] + translateVector[0]
        structureCoords[i][1] = structureCoords[i][1] + translateVector[1]
        structureCoords[i][2] = structureCoords[i][2] + translateVector[2]
    
# Calculate the shortest distance between two structures
def shortestDistance( peptideCoord, ligandCoord ):
    shortestDist = 100
    peptideAtomID = 0
    ligandAtomID = 0
    curDistance = 0
    for i in peptideCoord:
        for n in ligandCoord:
            curDistance = distance( peptideCoord[i], ligandCoord[n] )
            if curDistance < shortestDist:
                shortestDist = curDistance
                peptideAtomID = i
                ligandAtomID = n
    shortestVector = calcVector( peptideCoord[ peptideAtomID ], ligandCoord[ ligandAtomID ] )
    return shortestVector

# Shorten the magnitude of a vector. Pass by reference    
def reduceVectorMagnitude( vector , magReduction ):
    mag = sum(x**2 for x in vector) ** .5
    reductionFactor = (mag-magReduction) / mag
    newVector = [x * reductionFactor for x in vector] 
    return newVector

# Calculate ion position
def calcIonPosition( peptideCoM, peptideCoords, ligandCoM, ligandCoords ):
    
    # Translate ligand to be on top of peptide
    ligandToPeptideVector = calcVector( peptideCoM, ligandCoM )
    translateStructure( ligandCoords, ligandToPeptideVector )
    
    # Get random 3d vector with 100A length and translate ligand along it
    randomVector = make_rand_vector(3, 100) 
    translateStructure( ligandCoords, randomVector )
    
    # Get shortest distance between peptide and ion and translate along that
    # Stay a distance of 3.5AA away from the peptide though
    shortestDist = shortestDistance( peptideCoords, ligandCoords )
    shortestDist = reduceVectorMagnitude( shortestDist, 3.0 )
    translateStructure( ligandCoords, shortestDist )
    
    # Calculate the translation to be passed to the next step
    # i.e. the translation from starting to final position of the ligand
    requiredTranslation = calcVector( centerOfMass(ligandCoords), ligandCoM )    
    print "Required translationg of ligand in LEaP: ", requiredTranslation    
    
    # Return the ligand coords
    return requiredTranslation


