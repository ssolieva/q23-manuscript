import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import json
import pathlib
import os
path_to_models = '/Users/ssolieva/Desktop/comp_mount/colabfold/CH01_single_mutants/results/'
residues_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'] # mutations

overall_list_of_jsons = []
overall_list_of_pdbs  = []
for res in range(len(residues_list)):
    resi = residues_list[res]
    list_of_jsons = []
    list_of_pdbs  = []
    for seq in range(25):
        print(resi, seq)
        for model in range(1,6):
            for seed in range(10):
                for anyrank in range(51):
                    if anyrank <= 9:
                        json_path = f'{path_to_models}CH01_HC_{seq}{resi}_LC_scores_rank_00{anyrank}_alphafold2__multimer_v3_model_{model}_seed_00{seed}.json'
                        pdb_path  = f'{path_to_models}CH01_HC_{seq}{resi}_LC_unrelaxed_rank_00{anyrank}_alphafold2_multimer_v3_model_{model}_seed_00{seed}.pdb'
                        if os.path.isfile(json_path) == True:
                            print(pdb_path)
                            list_of_jsons.append(json_path)
                            list_of_pdbs.append(pdb_path)
                            #print('yes', json_path)
                    if anyrank > 9:
                        json_path = f'{path_to_models}CH01_HC_{seq}{resi}_LC_scores_rank_0{anyrank}_alphafold2_multimer_v3_model_{model}_seed_00{seed}.json'
                        pdb_path  = f'{path_to_models}CH01_HC_{seq}{resi}_LC_unrelaxed_rank_0{anyrank}_alphafold2_multimer_v3_model_{model}_seed_00{seed}.pdb'
                        if os.path.isfile(json_path) == True:
                            print(pdb_path)
                            list_of_jsons.append(json_path)
                            list_of_pdbs.append(pdb_path)
    overall_list_of_jsons.append(list_of_jsons)
    overall_list_of_pdbs.append(list_of_pdbs)
    
np.save("overall_list_of_jsons.npy", overall_list_of_jsons)
np.save("overall_list_of_pdbs.npy", overall_list_of_pdbs)
