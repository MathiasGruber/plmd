
DESIGN REQUIREMENTS
PLMD: Peptide Ligand Molecular Dynamics
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
* Upon completion be able to send plots & pdb files in zip package to email specified in config file
* Should contain a module for visualizing the molecules, and create images & vmd states of all the different stages automatically
