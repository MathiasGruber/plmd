PLMD: Peptide Ligand Molecular Dynamics
==============

The purpose of this module is automatize the use of the Molecular Dynamics 
framework Amber12 for investigation of the interaction between specifically 
small peptides and their ligands. Given the peptide and its ligands as input 
files, the module aim to make the process of setting up & running the molecular 
dynamics simulations and the subsequent analysis as easy as possible.

The modules furthermore aims to be easily extensible, so that as more and more 
simulation setups needs to be explored, these can easily be integrated into the
PLMD module. 

PLMD REQUIREMENTS
====================

* PLMD requires both Amber12 and AmberTools to be installed and set up properly in 
environment.

* The scripts were specifically designed to run at the DTU HPC cluster: http://www.cc.dtu.dk, but may work elsewhere, with appropriate modifications, as well


Installation
============================

You can get the latest version of this git repository by following command:

```
git clone git://github.com/MathiasGruber/plmd.git
```

This will create a "plmd" directory with all the required files. Before you can use the package it
must be linked within your environment. This can be done by running following commands in your shell:

```
PLMDHOME=path-to-plmd
export PLMDHOME
PATH=$PLMDHOME/bin:$PATH
export PYTHONPATH=$PYTHONPATH:$PLMDHOME/src
export PATH
```

where path-to-plmd represents the location of the plmd/ folder. This 
snippet can conveniently be placed at the bottom of your .bashrc file, so that
it'll be available on all future sessions.

How To: Setting up a simulation
=============================

* Create a directory where you want to run your simulations
* Create/retrieve your input files (peptide & ligand)
* Create a configuration file in your directory (see plmd/tools/defaultExample.cfg for sample)
* Run the command "plmd_setup configFile.cfg" to setup case directories
* Type "plmd_setup -h" for further information

How To: Submit to queue
=============================================

* To submit a case directory to the cluster queue use: "plmd_submit cases/" where cases/ is the directory created during previous step.
* The plmd_submit command will walk through all subdirectories and submit all identified cases.

How To: Run Post-processing
===========================

* To run postprocessing on a given case use the command "plmd_analyze cases/" which will run postprocessing on all cases found in your cases/ directory.
* Type "plmd_analyze -h" for further information

