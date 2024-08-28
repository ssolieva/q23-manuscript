import numpy as np
import os.path
import argparse, os
from frequency_toolbox import *

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
print(mount_path)

save_path = f'{mount_path}data/parsed_logfile/'
load_path = save_path

cdrh3_seqs_by_donor       = np.load(f'{load_path}cdrh3_seqs_by_donor.npy', allow_pickle=True)
donors_list = np.load(f'{load_path}donors_list.npy')

frequency_count = generate_frequency_before_removing_redundant_seqs(donors_list, cdrh3_seqs_by_donor, save_path)


donors_over1M = ['COS4', 'Donor-1', 'Donor-HD1', 'Donor-HD3', 'HD1', 'HD2', 'Subject-316188', 'Subject-326650', 'Subject-326651', 'Subject-326713', 'Subject-326737', 'Subject-326780', 'Subject-326797', 'Subject-326907', 'Subject-327059', 'Subject-BD1+3+4', 'Subject-BD1', 'Subject-BD14', 'Subject-BD1h', 'Subject-BD3', 'Subject-BD5', 'Subject-BD5h', 'Subject-BD6', 'Subject-CAP301', 'Subject-CAP312', 'Subject-CAP322', 'Subject-CAP351', 'Subject-CB2', 'Subject-CORD1', 'Subject-CORD2', 'Subject-CORD3', 'Subject-D103', 'Subject-D149', 'Subject-D181', 'Subject-HIP1', 'Subject-HIP2', 'Subject-HIP3']
donors_to_check = np.unique(donors_list)
list_donors_1M = []
for i in range(len(donors_to_check)):
    if donors_to_check[i] in donors_over1M:
        #print(i)
        list_donors_1M.append(donors_to_check[i])

print(len(list_donors_1M))


