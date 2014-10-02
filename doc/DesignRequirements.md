
== DESIGN REQUIREMENTS
== PLMD: Peptide Ligand Molecular Dynamics
==========================================

==== 1. Simulation Setup Module
===============================

* Easy plug & play of multiple peptides & ions into MD simulations using Amber12 and AmberTools
* Set up multiple cases with random ligand posisions
* Easy Control the details of the Amber input files via config file
* Easy Control the details of the server submission scripts via config file
* Save appropriate files required for subsequent analysis

==== 2. Simulation Result Analysis
==================================

* Automatically merge related MD trajectories from MD runs using `cpptraj` utility
* Convert trajectories from .mdcrd into binpos format
* Create plots of all common MD information; 
RMSd, energies, dihedral angles, end-to-end distance, radius of gyration, b factor analysis,
