import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pandas as pd
from yes_toolbox_turntypeIIIincluded import *

pdb_paths_and_cdrh3_seqs = []
with open(f'pdb_paths_with_cdrh3_seqs_new.txt') as f:
    pdb_paths_and_cdrh3_seqs.append([line.rstrip() for line in f])
    
pdb_paths = []
cdrh3_seqs = []
for i in range(len(pdb_paths_and_cdrh3_seqs[0])):
    pdb_paths.append(pdb_paths_and_cdrh3_seqs[0][i].split(',')[0])
    cdrh3_seqs.append(pdb_paths_and_cdrh3_seqs[0][i].split(',')[1])

beta_sheet_distance_cutoff = 20 #20A
for i in range(len(cdrh3_seqs)):
    pdb_path = pdb_paths[i]
    name = pdb_path[119:-4]
    print(i, name)
    cdrh3_sequence = cdrh3_seqs[i]
    phi_angles, psi_angles, dssp_values, hbonds3, hbonds4, distances = run_all_commands(name, pdb_path, cdrh3_sequence)
    features = [cdrh3_sequence, phi_angles, psi_angles, dssp_values, hbonds3, hbonds4, distances]
    csv_file = generate_dataframe(features, name) #  saving csv .
    turns_index,   beta_index    = classify_and_combine(csv_file)
    hbonds = classify_hbonds(csv_file['hbondi_i3'], csv_file['hbondi_i4'])
    combined_TB = combine_all_dist(beta_index, turns_index, hbonds, distances,beta_sheet_distance_cutoff)
    classification_of_hh = classify_turn_beta_distance_v2(csv_file, name, combined_TB, beta_sheet_distance_cutoff)
    if classification_of_hh == 1:
        print(f'yesHH,{pdb_path},{cdrh3_sequence}')
        plot_turn_beta_distance(csv_file, name, combined_TB, beta_sheet_distance_cutoff)
    if classification_of_hh == 0:
        print(f'noHH,{pdb_path},{cdrh3_sequence}')
        plot_turn_beta_distance(csv_file, name, combined_TB, beta_sheet_distance_cutoff)

