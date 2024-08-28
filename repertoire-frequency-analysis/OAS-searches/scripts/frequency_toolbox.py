def generate_frequency_before_removing_redundant_seqs(donors_list, cdrh3_seqs_by_donor, save_path):
    '''
    All data that goes here must be before redundant seqs are removed
        cdrh3_seqs_by_donor = no seqs were removed. 
    '''
    import numpy as np
    donors_list_unique = np.unique(donors_list)
    print('calculating frequencies -- before any redundant sequences are removed')
    frequency_count = []
    for donors in range(len(donors_list_unique)):
        donor = donors_list_unique[donors]
        n_seq_hits = len(cdrh3_seqs_by_donor[donors]) 
        
        # get the number of sequences from each donor in the full OAS database
        counter = []
        with open(f"/Users/ssolieva/Desktop/Kulp_lab/projects/OAS_database_searches/donor_counts/{donor}.out", "r") as fd: # import the results file
            file_contents = fd.read().splitlines() # read in the results
            for i in range(len(file_contents)):
                counter.append(int(file_contents[i]))
            #print(f'{donor}, Number of sequences in OAS database: {sum(counter)}')
            n_seqs_donor = sum(counter)
        fd.close()
        
        # calculate frequency per million for each donor
        if n_seqs_donor > 0:
            frequency_count.append(1000000*(n_seq_hits/n_seqs_donor))
        if n_seqs_donor == 0:
            print(n_seqs_donor)
            frequency_count.append(0) # remember: 0 will not show up on a log plot 
    np.save(f"{save_path}frequency_related/frequency_counts_before_removing_redundant_sequences.npy", frequency_count)
    return frequency_count 
