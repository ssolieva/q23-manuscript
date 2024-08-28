import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pandas as pd
from HH_toolbox_turntypeIIIincluded import *

save_path = '/export/home/WG-shahlo/projects/HIV/simulations/simulation_analysis/HH-ratios/data_every10/'

import argparse, os
# add argument for path
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser(prog='classification',
                                 description='classify frames as axe or not axe')
parser.add_argument('seqName',
                    type=str,
                    help=' ')
parser.add_argument('simPath',
                    type=dir_path,
                    help='Path to directory containing the simulation runs')
parser.add_argument('cdrh3Seq',
                    type=str,
                    help=' ')
parser.add_argument('run',
                    type=str,
                    help=' ')
args = parser.parse_args()

seq_name = args.seqName
sim_dir = args.simPath
cdrh3_sequence = args.cdrh3Seq
which_run = args.run

print(seq_name, sim_dir, cdrh3_sequence, which_run)

pdb_path = f'{sim_dir}input_files/{seq_name}-prot-masses.pdb'
traj = md.load(f'{sim_dir}/{which_run}/frame0_masses.xtc', top=pdb_path)

beta_sheet_distance_cutoff = 20 #20A
classification_list=[]
for index_n in np.arange(0,len(traj),10):##range(len(traj)):
    structure = md.load_frame(f'{sim_dir}/{which_run}/frame0_masses.xtc', index=index_n, top=pdb_path)
    name = f'{seq_name}'
    phi_angles, psi_angles, dssp_values, hbonds3, hbonds4, distances = run_all_commands(name, pdb_path, cdrh3_sequence, structure)
    features = [cdrh3_sequence, phi_angles, psi_angles, dssp_values, hbonds3, hbonds4, distances]
    csv_file = generate_dataframe(features) # not saving csv for space.
    turns_index,   beta_index    = classify_and_combine(csv_file)
    hbonds = classify_hbonds(csv_file['hbondi_i3'], csv_file['hbondi_i4'])
    combined_TB = combine_all_dist(beta_index, turns_index, hbonds, distances,beta_sheet_distance_cutoff)
    classification_of_hh = classify_turn_beta_distance(csv_file, name, combined_TB, beta_sheet_distance_cutoff)
    #nametitle = f'{name}{which_run}_frame{index_n}'
    #plot_turn_beta_distance(csv_file, nametitle, combined_TB, beta_sheet_distance_cutoff)
    #print(index_n, features)
    if classification_of_hh == 1:
        classification_list.append(1)
    if classification_of_hh == 0:
        classification_list.append(0)
np.save(f'{save_path}every10frames_classification_list_{seq_name}_{which_run}.npy', classification_list)

print(sum(classification_list), len(classification_list))
print('ratio:',sum(classification_list)/len(classification_list))
