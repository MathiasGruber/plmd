PLMD DESIGN REQUIREMENTS
==========================================

The two major components of this package are a module for setting up & running multiple molecular dynamics simulations of a given peptide-ligand interaction, and a module for subsequent analysis of the results.

1. Simulation Setup Module
===============================

* Easy plug & play of multiple peptides & ions into MD simulations using Amber12 and AmberTools
* Set up multiple cases with random ligand posisions
* Easy Control the details of the Amber input files via config file
* Easy Control the details of the server submission scripts via config file
* Save appropriate files required for subsequent analysis

2. Simulation Result Analysis
==================================

* Automatically merge related MD trajectories from MD runs using `cpptraj` utility
* Convert trajectories from .mdcrd into binpos format
* Create plots of following MD information; 
  1. RMSd (on different reference points), 
  2. Potential, kinetic and total energy, 
  3. Psi, phi and omega angles, 
  4. End-to-end distance, 
  5. Radius of gyration, 
  6. B factor analysis,
  7. Ligand to residue distance plots
  8. Time correlation
  9. 2d RMSd analysis
  10. Cluster analysis of peptide structures
  11. Hydrogen binding
  12. Block analysis
* Upon completion be able to send plots & pdb files in zip package to ftp specified in config file
* Should contain a module for visualizing the molecules, and create images & vmd states of all the different stages automatically

## 2.2. Design of Analysis Module
The main focus of this package is on the analysis module, which is expected to be continously updated with new features & improvements. The design of this module therefore warrents further explanation. 

The analysis modules are distributed in 3 directories in `src/plmd` called `caseAnalyses`, `globalAnalyses`, and `generalAnalyses`. As their names indicate, the first two refer to the analysis of single cases, and the comparison of multiple cases, respectively. `generalAnalyses` contain modules used by analyses in both single and multi-case analyses. 

The main class that manages the analyses is the one set up in `caseAnalysis.py`. This class manages which type (single or global) analysis the user wants to do, and it handles common tasks such as merging of trajectories, and passes all the neccesary information to the `__init__.py` file in either `caseAnalyses` or `globalAnalyses`, which then handles the actual analysis.

The package makes heavy use of the data processing power of cpptraj (supplied with Amber). These are mainly run by the `__init__.py` in `caseAnalyses`, in three runs, using the `cpptraj` templates `full`, `short`, `global` which are found in `plmd/templates` folder. Anyone wishing to add additional cpptraj processing should make use of the appropriate template, so as to minimize the amount of individual cpptraj runs.
