### This repository contains code from our manuscript, "Deep Mining of the Human Antibody Repertoire Identifies Frequent and Immunogenetically Diverse CDRH3 Topologies Targetable by Vaccination"

## Outline:
	-> MD-simulations
	-> insilico-saturated-mutagenesis
	-> pymol_session
	-> repertoire-frequency-analysis

### MD-simulations
	-> analysis
	  -> data
	  -> scripts
	    -> axe_classification_every10frames.py
	    -> axe_toolbox_turntypeIIIincluded.py
	    -> submit_slurm.sh
	  -> sim-analysis.ipynb
	-> input-files
	  -> Makefile
	  -> charmm-250ns.mdp
	  -> ions.mdp
	  -> minim.mdp
	  -> npt-charmm.mdp
	  -> nvt-charmm.mdp
	-> starting-structures
	  -> 3tcl_ch04.pdb
	  -> Ab41328.pdb
          -> V033mat.pdb
          -> seq10821.pdb
          -> seq10941.pdb
          -> seq12552.pdb
          -> seq15061.pdb
          -> seq17216.pdb
          -> seq20621.pdb
          -> seq215.pdb
          -> seq22299.pdb

### insilico-saturated-mutagenesis
	-> calculating_and_plotting_rsmd_ch01_mature_control_vs_3tcl.ipynb
	-> saturated-mutagenesis-analysis.ipynb
	-> save_paths_to_models.py

### pymol_sessions
	-> Figure_1
	-> Figure_2
	-> Figure_S2

### repertoire-frequency-analysis
	-> OAS-searches
	-> alphafold2-prediction-files/scripts
	-> axe-structure-analysis


