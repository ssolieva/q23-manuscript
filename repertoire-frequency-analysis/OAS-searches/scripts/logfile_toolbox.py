def open_all_log(logfile_path):
    '''
    open and read the log file and print some preview info if preview=True. 
    '''
    # open the file:
    with open(logfile_path, "r") as fd: # import the results file 
        logfile = fd.read().splitlines() # read in the results
        logfile_split = [] 
        for line in logfile:
            logfile_split.append(line.split()) # save each line to the list logfile_split
    fd.close()
    return logfile_split


def save_names_sequences_genes(logfile_split):
    donor_names = []
    cdrh3_seqs  = []
    full_seqs   = []
    dgene = []
    jgene = []
    vgene = []
    for line in range(len(logfile_split)): # for each line 
        for entry in range(len(logfile_split[line])): # for each entry in each line
            if logfile_split[line][entry].rfind("IGHV") == 0: # find which entry the Vgene is at 
                vgene.append(logfile_split[line][entry]) # save the v gene
                dgene.append(logfile_split[line][entry+1]) # save the d gene 
                jgene.append(logfile_split[line][entry+2]) # save the j gene
                if 'csv' in logfile_split[line][entry-1]: # if a .csv file is right before the Vgene, 
                    donor_names.append("no_donor_name_found") # then this line has no donor name
                else: # if the .csv line is not right before the Vgene, 
                    donor_names.append(logfile_split[line][entry-1]) # then save the reported donor name  
            if logfile_split[line][entry].isdigit() == True:
                if int(logfile_split[line][entry]) == len(logfile_split[line][entry-1]):
                    cdrh3_seqs.append(logfile_split[line][entry-1])     
        full_seqs.append(logfile_split[line][-2])
    return donor_names, cdrh3_seqs, full_seqs, dgene, jgene, vgene

def remove_no_and_no_name(donor_names, cdrh3_seqs, full_seqs, dgene, jgene, vgene, save_path):
    import numpy as np
    updated_donor_names = []
    updated_cdrh3_seqs = []
    updated_dgene = []
    updated_full_seqs = [] 
    updated_vgene = []
    updated_jgene = []
    for i in range(len(donor_names)):
        if donor_names[i] != 'no_donor_name_found' and donor_names[i] != 'no':
            updated_donor_names.append(donor_names[i])
            updated_cdrh3_seqs.append(cdrh3_seqs[i])
            updated_dgene.append(dgene[i])
            updated_full_seqs.append(full_seqs[i])
            updated_jgene.append(jgene[i])
            updated_vgene.append(vgene[i])
    np.save(f'{save_path}updated_donor_names', updated_donor_names)
    np.save(f'{save_path}updated_cdrh3_seqs', updated_cdrh3_seqs)
    np.save(f'{save_path}updated_full_seqs', updated_full_seqs)
    np.save(f'{save_path}updated_dgene', updated_dgene)
    np.save(f'{save_path}updated_jgene', updated_jgene)
    np.save(f'{save_path}updated_vgene', updated_vgene)
    return updated_donor_names, updated_cdrh3_seqs, updated_full_seqs, updated_dgene, updated_jgene, updated_vgene

def group_according_to_donors(updated_donor_names, updated_cdrh3_seqs, updated_full_seqs, save_path):
    import numpy as np
    full_seqs = updated_full_seqs
    cdrh3_seqs = updated_cdrh3_seqs
    donor_names = updated_donor_names
    unique_donor_names = np.unique(updated_donor_names)

    cdrh3_seqs_by_donor = []
    full_seqs_by_donor = []
    for i in range(len(unique_donor_names)):
        tmp_cdrh3 = []
        tmp_fullseq = []
        for j in range(len(donor_names)):
            tmp_tmp_cdrh3 = []
            tmp_tmp_fullseq = []
            if donor_names[j] == unique_donor_names[i] and donor_names[j] != 'no' and donor_names[j] != 'no_donor_name_found':
                tmp_tmp_cdrh3.append(cdrh3_seqs[j])
                tmp_tmp_fullseq.append(full_seqs[j])
            tmp_cdrh3.append(tmp_tmp_cdrh3)
            tmp_fullseq.append(tmp_tmp_fullseq)
        cdrh3_seqs_by_donor.append(np.concatenate(tmp_cdrh3))
        full_seqs_by_donor.append(np.concatenate(tmp_fullseq))

    cdrh3_seqs_by_donor_arr = np.array(cdrh3_seqs_by_donor, dtype=object)
    full_seqs_by_donor_arr = np.array(full_seqs_by_donor,dtype=object)
    np.save(f'{save_path}cdrh3_seqs_by_donor',cdrh3_seqs_by_donor_arr)
    np.save(f'{save_path}full_seqs_by_donor',full_seqs_by_donor_arr)
    return unique_donor_names, cdrh3_seqs_by_donor, full_seqs_by_donor

def generate_unique_seqs(cdrh3_seqs_by_donor, full_seqs_by_donor, save_path):
    import numpy as np
    unique_cdrh3_seqs_by_donor_ind = []
    unique_cdrh3_seqs_by_donor = []
    for i in range(len(cdrh3_seqs_by_donor)):
        unique_cdrh3_seqs_by_donor.append(np.unique(cdrh3_seqs_by_donor[i], return_index=True)[0])
        unique_cdrh3_seqs_by_donor_ind.append(np.unique(cdrh3_seqs_by_donor[i], return_index=True))
        
    unique_full_seq_by_donor = []
    for i in range(len(full_seqs_by_donor)):
        unique_full_seq_by_donor_temp = []
        for j in range(len(unique_cdrh3_seqs_by_donor_ind[i][1])):
            n = unique_cdrh3_seqs_by_donor_ind[i][1][j]
            unique_full_seq_by_donor_temp.append(full_seqs_by_donor[i][n])
        unique_full_seq_by_donor.append(unique_full_seq_by_donor_temp)

    unique_cdrh3_seqs_by_donor_arr = np.array(unique_cdrh3_seqs_by_donor, dtype=object)
    unique_full_seq_by_donor_arr = np.array(unique_full_seq_by_donor, dtype=object)
    np.save(f'{save_path}unique_cdrh3_seqs_by_donor', unique_cdrh3_seqs_by_donor_arr)
    np.save(f'{save_path}unique_full_seq_by_donor'  , unique_full_seq_by_donor_arr)
    return unique_cdrh3_seqs_by_donor, unique_full_seq_by_donor

def diff_letters(a,b):
    '''reports the number of different letters (residues)'''
    return sum ( a[i] != b[i] for i in range(len(a)) )

def remove_almost_duplicates(unique_cdrh3_seqs_by_donor, unique_full_seq_by_donor, save_path):
    '''
    data = unique_cdrh3_seqs_by_donor[0]
    '''
    import numpy as np
    dataset = unique_cdrh3_seqs_by_donor
    dataset_full_seq = unique_full_seq_by_donor 
    
    different_seqs_final      = [] # final cdrh3 sequences
    different_seqs_final_full = [] # final full sequences
    
    # get unique cdrh3 lengths:
    unique_lengths = []
    for seq in range(len(dataset)):
        unique_lengths.append(len(dataset[seq]))
    unique_lengths = np.unique(unique_lengths)
    
    # group data into unique length groups:
    unique_length_groups = []
    unique_length_groups_full = []
    for length in range(len(unique_lengths)):
        unique_length_groups_0 = []
        unique_length_groups_0_full = []
        for seq in range(len(dataset)):
            unique_length_groups_00 = []
            unique_length_groups_00_full = []
            if len(dataset[seq]) == unique_lengths[length]:
                unique_length_groups_00.append(dataset[seq])
                unique_length_groups_00_full.append(dataset_full_seq[seq])
            unique_length_groups_0.append(unique_length_groups_00)
            unique_length_groups_0_full.append(unique_length_groups_00_full)
        unique_length_groups.append(np.concatenate(unique_length_groups_0))
        unique_length_groups_full.append(np.concatenate(unique_length_groups_0_full))
    #print(d, len(unique_length_groups), unique_lengths)
    
    # for each unique length group, get rid of duplicates and almost duplicates
    for group in range(len(unique_length_groups)):
        if len(unique_length_groups[group]) == 1: # if there's only 1 sequence, then save that sequence
            #print(len(unique_length_groups[group]), unique_length_groups[group][0])
            different_seqs_final.append(unique_length_groups[group][0])
            different_seqs_final_full.append(unique_length_groups_full[group][0])
        
        all_seqs_updated = []
        all_seqs_updated_full = [] # new line
        
        if len(unique_length_groups[group]) > 1: # if more than 1 sequence, 
            all_seqs = unique_length_groups[group]
            all_seqs_full = unique_length_groups_full[group] #unique_length_groups[group] # new line
            
            a = unique_length_groups[group]
            b = unique_length_groups[group]
            pairs = [(i, j) for i in a for j in b if i < j] # make pairwise pairs of the sequences
            #print('number of pairs:', len(pairs))
            
            for pair in range(len(pairs)):
                n_diff = diff_letters(pairs[pair][0], pairs[pair][1])
                if n_diff <= (int(len(pairs[pair][0])*0.15)):                 # 15% / 85%
                   # print(n_diff, pairs[pair][0], pairs[pair][1])
                    result = np.where(all_seqs==f'{pairs[pair][1]}')[0]
                    all_seqs = np.delete(all_seqs, result) 
                    all_seqs_full =  np.delete(all_seqs_full, result) 
            for al in range(len(all_seqs)):
                all_seqs_updated.append(all_seqs[al])
                all_seqs_updated_full.append(all_seqs_full[al])
        #print('all_seqs_updated', all_seqs_updated[0],  all_seqs_updated_full[0])
        #print(all_seqs_updated, len(all_seqs_updated))
        if len(all_seqs_updated) != 0:
            for ds in range(len(all_seqs_updated)):
                different_seqs_final.append(all_seqs_updated[ds])
                different_seqs_final_full.append(all_seqs_updated_full[ds])
    redundant_seqs_removed_within_donors_cdrh3 = different_seqs_final
    redundant_seqs_removed_within_donors_fullseq = different_seqs_final_full
    return redundant_seqs_removed_within_donors_cdrh3, redundant_seqs_removed_within_donors_fullseq








