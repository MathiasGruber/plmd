plmd
====

PLMD: Peptide Ligand Molecular Dynamics, is a python module for managing the creation, running and analysis of molecular dynamics simulations of peptide-ion interactions using Amber12 and AmberTools.

== DESCRIPTION 
== PLMD: Peptide Ligand Molecular Dynamics
==============

The purpose of this module is automatize the use of the Molecular Dynamics 
framework Amber12 for investigation of the interaction between specifically 
small peptides and their ligands. Given the peptide and its ligands as input 
files, the module aim to make the process of setting up & running the molecular 
dynamics simulations and the subsequent analysis as easy as possible.

The modules furthermore aims to be easily extensible, so that as more and more 
simulation setups needs to be explored, these can easily be integrated into the
PLMD module. 

== Installation
============================

To install the library, you must save the PLMD/ folder to a location of your 
choice, and then specify the following environment variable in your shell:

PLMDHOME=/zhome/20/f/54098/PLMD
export PLMDHOME
PATH=$PLMDHOME/bin:$PATH
export PYTHONPATH=$PYTHONPATH:$PLMDHOME/src
export PATH

where /zhome/20/f/54098/PLMD represents the location of the PLMD folder. This 
snippet can conveniently be placed at the bottom of your .bashrc file, so that
it'll be available on all future bash sessions.

== PLMD REQUIREMENTS
====================

PLMD requires the module "mdanalysis" to work. This module can be retrieved at:
https://code.google.com/p/mdanalysis/


== Setting up a simulation
=============================

* Create a directory where you want to run your simulations
* Create/retrieve your input files (peptide & ligand)
* Run the command "setupAmberLigandSim peptide.pdb ligand.mol2"
* This will create the file setupAll.py in your directory. This file contains all the information required to setup the simulation

== Phase 2: Setup the cases & submit to queue
=============================================

* Once you're done editing the setupAll.py file, you execute it by typing "python setupAll.py"
* This will set up all the specified cases, input files etc. 
* It will also create a file called "submitAll.py".
* Running "python submitAll.py" will submit all the cases to the queue.

== Phase 3: Post-processing
===========================

* Once the simulations are completed

