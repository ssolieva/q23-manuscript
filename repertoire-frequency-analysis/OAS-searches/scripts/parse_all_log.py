import numpy as npl
import matplotlib.pyplot as plt
import os.path
import math
import argparse, os
from logfile_toolbox import *

# add arguments
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
parser = argparse.ArgumentParser(prog='parse-logfile',
                                 description='parse through a logfile.txt from OAS database searches.')
parser.add_argument('path', 
                    type=dir_path, 
                    help='Path to directory containing the logfile.txt')
args = parser.parse_args()

mount_path = args.path

import numpy as np
import matplotlib.pyplot as plt
import os.path
import math
import argparse, os
from logfile_toolbox import *

save_path = f'{mount_path}data/parsed_logfile/'
load_path = save_path
logfile_path = f'{mount_path}logfile/all_log.txt'

logfile_split = open_all_log(logfile_path)
donor_names, cdrh3_seqs, full_seqs, dgene, jgene, vgene = save_names_sequences_genes(logfile_split)
updated_donor_names, updated_cdrh3_seqs, updated_full_seqs, updated_dgene, updated_jgene, updated_vgene = remove_no_and_no_name(donor_names, cdrh3_seqs, full_seqs, dgene, jgene, vgene, save_path)
unique_donor_names, cdrh3_seqs_by_donor, full_seqs_by_donor = group_according_to_donors(updated_donor_names, updated_cdrh3_seqs, updated_full_seqs, save_path)
unique_cdrh3_seqs_by_donor, unique_full_seq_by_donor = generate_unique_seqs(cdrh3_seqs_by_donor, full_seqs_by_donor, save_path)

redundant_seqs_removed_within_donors_cdrh3_final = []
redundant_seqs_removed_within_donors_fullseq_final = []
for donor in range(len(unique_cdrh3_seqs_by_donor)):
    redundant_seqs_removed_within_donors_cdrh3, redundant_seqs_removed_within_donors_fullseq = remove_almost_duplicates(unique_cdrh3_seqs_by_donor[donor], unique_full_seq_by_donor[donor], save_path)
    redundant_seqs_removed_within_donors_cdrh3_final.append(redundant_seqs_removed_within_donors_cdrh3)
    redundant_seqs_removed_within_donors_fullseq_final.append(redundant_seqs_removed_within_donors_fullseq)

np.save(f'{save_path}redundant_seqs_removed_within_donors_cdrh3_final', redundant_seqs_removed_within_donors_cdrh3_final)
np.save(f'{save_path}redundant_seqs_removed_within_donors_fullseq_final', redundant_seqs_removed_within_donors_fullseq_final)

np.save(f'{save_path}redundant_seqs_removed_within_donors_cdrh3_final_concat'  , np.concatenate(redundant_seqs_removed_within_donors_cdrh3_final))
np.save(f'{save_path}redundant_seqs_removed_within_donors_fullseq_final_concat', np.concatenate(redundant_seqs_removed_within_donors_fullseq_final))

updated_donor_names       = np.load(f'{load_path}updated_donor_names.npy')
unique_donor_names = np.unique(updated_donor_names)
donors_list = []
for i in range(len(redundant_seqs_removed_within_donors_cdrh3_final)):
    for j in range(len(redundant_seqs_removed_within_donors_cdrh3_final[i])):
        donors_list.append(unique_donor_names[i])
np.save(f'{save_path}donors_list.npy', donors_list) 

# print statements:
print('----------------------------------------------------------------------------------')
print(mount_path)
print('\t>Number of donors:',len(unique_donor_names))
print('\t>Number of sequences before redundancy removed:', len(cdrh3_seqs))
print('\t>Number of sequences after redundancy removed:',len(np.concatenate(redundant_seqs_removed_within_donors_cdrh3_final)))
print('----------------------------------------------------------------------------------')
