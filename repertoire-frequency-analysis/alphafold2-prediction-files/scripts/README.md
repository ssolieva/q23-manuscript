# Scripts to set up and run AlphaFold2 predictions.

        ---> example-convert_npy_to_a3m.py

This is an example script for how numpy files containing variable HC sequences were converted into a3m files. a3m files are the custom t
emplates given for AlphaFold2 predictions. 

	---> pre_gpu_custom_template.sh

This script submits AlphaFold2 predictions to the cluster. Use the following command for submission:
for d in dir*; do cd $d; sbatch ../pre_gpu_custom_template.sh; cd -; done;

	---> cass_path 
		---> cass.cif 
		---> pdb70_a3m.ffdata
		---> pdb70_a3m.ffindex
		---> pdb70_cs219.ffdata
		---> pdb70_cs219.ffindex

This directory is called by pre_gpu_custom_template.sh, it contains information needed for the AlphaFold2 predictions. cass.cif is an antibody structure with its CDRH3 removed. This antibody structure is CH01, PG9, or V033 depending on which set of sequences was being predicted.
