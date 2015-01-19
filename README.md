PLMD: Peptide Ligand Molecular Dynamics
==============

The purpose of this module is automatize the use of the Molecular Dynamics 
framework Amber14 for investigation of the interaction between small peptides 
and their ligands. Given the peptide and its ligands as input 
files, the module aim to make the process of setting up & running the molecular 
dynamics simulations and the subsequent analysis as easy as possible.

The modules furthermore aims to be easily extensible, so that as more and more 
simulation setups needs to be explored, these can easily be integrated into the
PLMD module. 

PLMD REQUIREMENTS
====================

* PLMD requires both Amber14 and AmberTools14 to be installed and set up properly in 
environment.

* The analysis module requires `MDAnalysis` to be installed:
```
https://code.google.com/p/mdanalysis/
```
See below on how to install this locally on the DTU HPC

* The scripts were specifically designed to run at the DTU HPC cluster: http://www.cc.dtu.dk, but may work elsewhere, with appropriate modifications, as well

* The module furthermore requires numpy, matplotlib, and pylab modules to be properly installed. As such, on the DTU HPC, it requires that you are logged in on the backend nodes!

* The module requires python 2.7. On the DTU HPC, python 2.7 can be loaded by typing the following:
```
module load python/2.7.3
```
This line can be added to your `.bashrc` file, so that you won't have to run it during every session. 

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

This module is still under constant development. Normal updates can be pulled with:

```
git pull
```

If this fails, e.g. in the case of local changes, you can update to the latest version (and delete all local changes) by doing a pull from within the plmd directory as follows:

```
git fetch --all
git reset --hard origin/master
``` 

## Installing MDAnalysis

The module requires MDAnalysis to be installed. This can be installed (on your local machine, or the HPC) by using the following commands:

```
git clone https://code.google.com/p/mdanalysis
cd mdanalysis
python setup.py build
python setup.py install --user
```
Please note that if you are installing it on the DTU HPC, you need to login on the backend, and use `module load python/2.7.3` first!

How To: Setting up a set of cases
=============================

* Create a directory where you want to run your simulations, e.g. "HPO4_SGAGKT"
* Go to that directory.
* Retrieve default configuration file with the command:
```
plmd_getConfig.py
```
* Review the configuration file. Ligand & peptide input files do not have to be specified.
* Setup cases according to config file by typing

```
plmd_setup.py defaultExample.cfg -pPeptide "NSER GLY ALA GLY LYS CTHR" -pLigand HPO4 --quiet
```
* This will create a cases/ directory with all your chases set up and ready for submission using the SGAGKT peptide, and the HPO4 ion. 
* Type "plmd_setup -h" for further information on setup options

How To: Submit to queue
=============================================

* To submit a case directory to the cluster queue use: 
```
plmd_submit defaultExample.cfg cases/
```
where `cases/` is the directory created during previous step and `defaultExample.cgf` is the configuration file
* The plmd_submit command will walk through all subdirectories and submit all identified cases. Submission specifications (walltime, nodes, job name etc.) are in the configuration file.

How To: Run Post-processing
===========================

* To postprocess a given case use the following command 
```
plmd_analyze defaultExample.cfg cases/
```
which will run postprocessing on all cases found in your cases/ directory.
* Type "plmd_analyze -h" for further information on the options available during the step. Some of the options warrent special notice:

`-noQueue:` This effectively runs the analysis in your terminal, instead of submitting it to the server queue.

`-noMerge` Each case is run in several separate MD runs. During analysis, these are merged together, which can take quite some time. If you don't need the latest merged trajectory (e.g. during debugging), you can use this option to turn merging off.

`-noBlock` The block averaging analysis inherently requires a lot of time. In case it's not strictly needed, you can disable it with this flag.

`-noFTP` By default the script will assume you want your results uploaded to a FTP server when they are done. The FTP details are specified in the configuration file (except password, which is not stored at any point). If you want to disable this feature, this flag is for you.

`-doGlobal` Once you've run your initial analysis, you may want to re-run the analysis module with this flag enabled. It will use data files from the previous analyses, and do a "global" analysis, i.e. compare various results from each case with each other.
 

How To: Add Analysis Module
===========================

If you are looking add a specific analysis module to this package, please read this first. There are two kinds of analysis modules: case-based and global analyses. The modules can be found in `src/plmd/caseAnalyses` and  `src/plmd/globalAnalyses`, respectively. When you want to add a new module, you'll have to add that in to these directories (see other modules for inspiration). The `__init__`-files in the module directories, are where you link your module to the rest of the package, i.e. make sure you:

* `import YOUR_PACKAGE` at the top of the `__init__`-file
* Add the call to your module in the `run_all`-function of the `__init__`-file

## Adding cpptraj analyses
Cpptraj is run through twice during the script; 

- First time only 500 frames from the trajectory are included (which is enough for plotting of data). If you want to add your analysis to this run, add it to `src/templates/cpptraj_analysis_short.txt`
- Second time the entire trajectory is read (which is suitable for H-bond analysis, cluster analysis, PCA and so forth. If you want to add your analysis to this run, add it to `src/templates/cpptraj_analysis_full.txt`

## Adding MDAnalysis analyses
In the `src/plmd/analyses/__init__`-file, in the constructor, the trajectory of the simulation is loaded and some basic MDAnalysis parameters are set up. These can be passed to your analysis module, so that this retrieval of data only has to happen once. See e.g. the `src/plmd/analyses/block.py` analysis module for how this is set up.
