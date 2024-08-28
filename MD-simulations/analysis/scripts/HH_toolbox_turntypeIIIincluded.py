import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from math import pi

# step 2: get CDRH3 index 
def generate_cdrh3_index(pdb_path, cdrh3_sequence):
    import mdtraj as md
    structure = md.load(pdb_path)
    all_ca = structure.topology.select('protein and name CA')
    rescode_list = []
    for i in range(len(all_ca)):
        rescode_list.append(structure.topology.atom(all_ca[i]).residue.code)
    rescodes = ''.join(rescode_list)
    cdrh3_index = rescodes.find(cdrh3_sequence)
    start_index = cdrh3_index
    end_index = cdrh3_index+len(cdrh3_sequence)
    if rescodes[start_index:end_index] != cdrh3_sequence:
        print(f'ERROR (generate_cdrh3_index) | actual seq: {cdrh3_sequence} | output seq: {rescodes[start_index:end_index]}')
    #print('\tsequence:',rescodes)
    #print('\tcdrh3 index:',cdrh3_index)
    #print(f'\tcdrh3 = {rescodes[start_index:end_index]}')
    return cdrh3_index

# step 3: get DSSP values
def generate_dssp(structure, cdrh3_sequence, cdrh3_index):
    '''
    generate dssp assignments for each residue given in the cdrh3_sequence input.
    '''
    import mdtraj as md
    dssp_orig = md.compute_dssp(structure, simplified=False)
    
    dssp = []
    for i in range(len(dssp_orig[0])):
        if dssp_orig[0][i] != ' ':
            dssp.append(dssp_orig[0][i])
        if dssp_orig[0][i] == ' ':
            dssp.append('L')
    
    #for i in range(len(dssp)):
    #    print(i, dssp[i], structure.topology.residue(i).code)
    
    start_index = cdrh3_index
    end_index = cdrh3_index+len(cdrh3_sequence)
    dssp_values = dssp[start_index:end_index]
    #print('DSSP:', dssp_values)
    
    return dssp_values

# step 4: get phi angles
def generate_cdrh3_phi_angles(structure, cdrh3_sequence, cdrh3_index):
    '''
    generate phi angles (degrees) for each residue given in the cdrh3_sequence input.
    phi_angles[1] = angles.
    phi_angles[0] = corresponding atoms.
    '''
    import mdtraj as md
    phi_angles = md.compute_phi(structure)
    start_index = cdrh3_index-1
    end_index = start_index + len(cdrh3_sequence)
    phi_angles_cdrh3 = phi_angles[1][0][start_index:end_index]

    #for i in range(len(phi_angles[0])):
    #    print(i, phi_angles[0][i][2])
    
    for i in range(len(cdrh3_sequence)):  # check that the cdrh3 residues are correct. 
        actual_residue = str(cdrh3_sequence[i])
        output_residue = str(structure.topology.atom(phi_angles[0][start_index+i][2]).residue.code)
        if actual_residue != output_residue:
            print(f'ERROR (generate_cdrh3_phi_angles) | actual residue = {actual_residue} | output residue = {output_residue}')
            
    return phi_angles_cdrh3

# step 5: get psi angles
def generate_cdrh3_psi_angles(structure, cdrh3_sequence, cdrh3_index):
    '''
    generate psi angles (degrees) for each residue given in the cdrh3_sequence input.
    psi_angles[1] = angles.
    psi_angles[0] = corresponding atoms.
    '''
    import mdtraj as md
    psi_angles = md.compute_psi(structure)
    start_index = cdrh3_index
    end_index = start_index + len(cdrh3_sequence)
    psi_angles_cdrh3 = psi_angles[1][0][start_index:end_index]
    #print(psi_angles[0][start_index])
    #for i in range(len(cdrh3_sequence)):
    #    output_residue = str(structure.topology.atom(psi_angles[0][start_index+i][1]).residue.code)
    #    print(cdrh3_sequence[i], psi_angles[0][start_index+i][1], output_residue)
    for i in range(len(cdrh3_sequence)):  # check that the cdrh3 residues are correct. 
        actual_residue = str(cdrh3_sequence[i])
        output_residue = str(structure.topology.atom(psi_angles[0][start_index+i][2]).residue.code)
        if actual_residue != output_residue:
            print(i,f'ERROR (generate_cdrh3_psi_angles) | actual residue = {actual_residue} | output residue = {output_residue}')
    return psi_angles_cdrh3

# step 6: check for i+3 and i+4 hydrogen bonds:
def generate_hbonds(structure, cdrh3_sequence, name, cdrh3_index):
    '''
    check for i+3 and i+4 hydrogen bond for each residue given in the cdrh3_sequence input.
    '''
    import mdtraj as md
    import numpy as np
    start_index = cdrh3_index
    end_index = cdrh3_index+len(cdrh3_sequence)
    all_inds = np.arange(start_index,end_index)
    N_atoms = []
    O_atoms = []
    for i in range(len(all_inds)):
        N_atom = structure.topology.select(f'resid {all_inds[i]} and name N')[0] # res index is 0-based.
        O_atom = structure.topology.select(f'resid {all_inds[i]} and name O')[0]
        N_atom_info = structure.topology.atom(N_atom)
        O_atom_info = structure.topology.atom(O_atom)
        N_atoms.append(N_atom)
        O_atoms.append(O_atom)
    hbond_i_to_i3 = []
    hbond3_resi_ind = []
    for i in range(len(N_atoms)-3):
        hbond_i_to_i3_opt1 = md.compute_distances(structure, [[N_atoms[i],  O_atoms[i+3]]])[0][0]*10 # nm to A
        hbond_i_to_i3_opt2 = md.compute_distances(structure, [[N_atoms[i+3],O_atoms[i]]])[0][0]*10 # nm to A
        #print(i, 'N to O:' ,hbond_i_to_i3_opt1, N_atoms[i], O_atoms[i+3])
        #print(i, 'O to N:',hbond_i_to_i3_opt2, N_atoms[i+3],O_atoms[i])
        if hbond_i_to_i3_opt1 <= 4.0 or hbond_i_to_i3_opt2 <= 4.0:
            #print('i to i+3 hbond:',structure.topology.atom(N_atoms[i]), structure.topology.atom(N_atoms[i+3]))
            hbond3_resi_ind.append(i)
        else:
            hbond3_resi_ind.append(0)
    for i in range(3): # to get full length to 20. 
        hbond3_resi_ind.append(0)
    hbond_i_to_i4 = []
    hbond4_resi_ind = []
    for i in range(len(N_atoms)-4):
        hbond_i_to_i4_opt1 = md.compute_distances(structure, [[N_atoms[i],  O_atoms[i+4]]])[0][0]*10 # nm to A
        hbond_i_to_i4_opt2 = md.compute_distances(structure, [[N_atoms[i+4],O_atoms[i]]])[0][0]*10 # nm to A
        if hbond_i_to_i4_opt1 <= 4.0 or hbond_i_to_i4_opt2 <= 4.0:
            #print('i to i+4 hbond:', structure.topology.atom(N_atoms[i]), structure.topology.atom(N_atoms[i+3]))
            hbond4_resi_ind.append(i)
        else:
            hbond4_resi_ind.append(0)
    for i in range(4):# to get full length to 20. 
        hbond4_resi_ind.append(0)
    return hbond3_resi_ind, hbond4_resi_ind
    

# step 7: generate distances from first resi to every resi in cdrh3

def generate_distances(structure, cdrh3_sequence, name, cdrh3_index):
    '''
    generate distances from first resi to every resi in cdrh3
    '''
    import mdtraj as md
    all_ca = structure.topology.select('protein and name CA')
    start_index = cdrh3_index
    end_index = cdrh3_index+len(cdrh3_sequence)
    cdrh3_ca_atoms = all_ca[start_index:end_index]
    #print(cdrh3_ca_atoms)
    for i in range(len(cdrh3_sequence)):
        #print(structure.topology.atom(cdrh3_ca_atoms[i]).residue.code)
        actual_residue = cdrh3_sequence[i]
        output_residue = structure.topology.atom(cdrh3_ca_atoms[i]).residue.code
        if actual_residue != output_residue:
            print(f'ERROR (generate_distances) | actual = {actual_residue} | output = {output_residue}')
    distance = []
    for i in range(len(cdrh3_ca_atoms)):
        distance.append(md.compute_distances(structure, [[cdrh3_ca_atoms[0], cdrh3_ca_atoms[i]]])[0][0]*10) # A
    return distance


# step 8: run all commands. 
def run_all_commands(name, pdb_path, cdrh3_sequence, structure):
    import mdtraj as md
    import numpy as np
    import matplotlib.pyplot as plt
    from math import pi
    
    cdrh3_index      = generate_cdrh3_index(pdb_path, cdrh3_sequence)
    dssp_values      = generate_dssp(structure, cdrh3_sequence, cdrh3_index)
    phi_angles       = generate_cdrh3_phi_angles(structure, cdrh3_sequence, cdrh3_index)*180/pi
    psi_angles       = generate_cdrh3_psi_angles(structure, cdrh3_sequence, cdrh3_index)*180/pi
    hbonds3, hbonds4 = generate_hbonds(structure, cdrh3_sequence, name, cdrh3_index)
    distances        = generate_distances(structure, cdrh3_sequence, name, cdrh3_index)
    
    return phi_angles, psi_angles, dssp_values, hbonds3, hbonds4, distances 

# Step 9: create a dataframe:

def generate_dataframe(features):
    ResidueIndex = []
    ResidueType  = []
    Phi = []
    Psi = []
    DSSP = []
    hbondi_i3 = []
    hbondi_i4 = []
    distances = []
    for i in range(len(features[0])):
        ResidueIndex.append(i)
        ResidueType.append(features[0][i])
        Phi.append(features[1][i])
        Psi.append(features[2][i])
        DSSP.append(features[3][i])
        hbondi_i3.append(features[4][i])
        hbondi_i4.append(features[5][i])
        distances.append(features[6][i])
    dataframe = {'ResidueIndex' : ResidueIndex, 'ResidueType' : ResidueType, 'Phi' : Phi, 'Psi' : Psi,
                 'DSSP' : DSSP, 'hbondi_i3': hbondi_i3, 'hbondi_i4': hbondi_i4, 'distances': distances}
    return dataframe





# classification / filtering the HH values:
def classify_phi_psi_TypeI(phi,psi):
    classify = []
    for i in range(len(phi)-2):
        if phi[i+1] > -100 and phi[i+1] < -25 and psi[i+1] > -75 and psi[i+1] < 25:
            if phi[i+2] > -150 and phi[i+2] < -50 and psi[i+2] > -50 and psi[i+2] < 50:
                classify.append(1)
            else:
                classify.append(0)
        else:
            classify.append(0)
    return classify

def classify_phi_psi_TypeII(phi,psi):
    # i+1 | Phi: -25 to -100 | Psi: 75 to 180
    # i+2 | Phi: 50 to 125 | Psi: -50 to 50
    classify = []
    for i in range(len(phi)-2): 
        if phi[i+1] > -100 and phi[i+1] < -25 and psi[i+1] > 75 and psi[i+1] < 180:
            if phi[i+2] > 50 and phi[i+2] < 125 and psi[i+2] > -50 and psi[i+2] < 50:
                classify.append(1)
            else:
                classify.append(0)
        else:
            classify.append(0)
    return classify

def classify_phi_psi_TypeIprime(phi,psi):
    # i+1 | Phi: 25 to 100 | Psi: 0 to 75
    # i+2 | Phi: 50 to 150 | Psi: -50 to 50
    classify = []
    for i in range(len(phi)-2): 
        if phi[i+1] > 25 and phi[i+1] < 100 and psi[i+1] > 0 and psi[i+1] < 75:
            if phi[i+2] > 50 and phi[i+2] < 150 and psi[i+2] > -50 and psi[i+2] < 60:
                classify.append(1)
            else:
                classify.append(0)
        else:
            classify.append(0)
    return classify

def classify_phi_psi_TypeIII(phi,psi):
    # i+1 | Phi: -25 to -100 | Psi: 75 to 180
    # i+2 | Phi: 50 to 125 | Psi: -50 to 50
    classify = []
    for i in range(len(phi)-2): 
        if phi[i+1] > -100 and phi[i+1] < -50 and psi[i+1] > -60 and psi[i+1] < 10:
            if phi[i+2] > -150 and phi[i+2] < -60 and psi[i+2] > 75 and psi[i+2] < 160:
                classify.append(1)
            else:
                classify.append(0)
        else:
            classify.append(0)
    return classify

def classify_phi_psi_beta(phi,psi):
    '''
    1 = beta sheet classification 
    0 = not classified as a beta sheet
    '''
    classify = []
    for i in range(len(phi)): # previously phi -65 
        if phi[i] < -50 and  psi[i] > 95:
            classify.append(1)
        else:
            if phi[i] < -50 and  psi[i] < -170:
                classify.append(1)
            else: 
                classify.append(0)
    return classify

def classify_and_combine(data_set):
    '''
    turns considers cdrh3 resis up until the last 2 resis.
    '''
    
    dssp_values = data_set['DSSP']
    
    typeI  = classify_phi_psi_TypeI(data_set['Phi'],      data_set['Psi'])
    typeII = classify_phi_psi_TypeII(data_set['Phi'],     data_set['Psi'])
    typeIp = classify_phi_psi_TypeIprime(data_set['Phi'], data_set['Psi'])
    typeIII = classify_phi_psi_TypeIII(data_set['Phi'], data_set['Psi'])
    typebeta = classify_phi_psi_beta(data_set['Phi'], data_set['Psi'])
   
    turns = []
    beta = []
    for i in range(len(typeI)):
        if dssp_values[i] in 'TS': #check dssp first
            turns.append(i)
        if dssp_values[i] not in 'TS':
            if typeI[i] == 1 or typeII[i] == 1 or typeIp[i] == 1 or typeIII[i] ==1: # including type 3 turn. 
                turns.append(i+1)
                turns.append(i+2)
    
    for i in range(len(dssp_values)):
        if dssp_values[i] == 'E' or dssp_values[i] == 'B':
            beta.append(i)
        if dssp_values[i] != 'E' and dssp_values[i] != 'B':
            if typebeta[i] == 1:
                beta.append(i)
    return turns, beta


def classify_hbonds(honds3, hbonds4):
    hbonds = []
    for i in range(len(honds3)):
        if honds3[i] != 0:
            hbonds.append(3)
        else:
            if hbonds4[i] !=0:
                hbonds.append(4)
            else: 
                hbonds.append(0)
    return hbonds

# check if it is a hammerhead or not according to the turn-beta-turn pattern:
def check_for_pattern(combined_data, distances, beta_sheet_distance_cutoff):
    '''Newest version, added on 07172024 '''
    combined_TB_joined = ''.join(combined_data)
    turn_resis = []
    beta_resis = []
    for i in range(len(combined_TB_joined)-2):
        if combined_TB_joined[i] in 'TFZ' and combined_TB_joined[i+1] in 'TFZ':
            turn_resis.append(i)
        if combined_TB_joined[i] in 'FD' and combined_TB_joined[i+1] in 'FD':
            beta_resis.append(i) # either far beta or far beta+turn

    # print(turn_resis, beta_resis)
    # non consecutive:
    from itertools import groupby
    def ranges(lst):
        pos = (j - i for i, j in enumerate(lst))
        t = 0
        for i, els in groupby(pos):
            l = len(list(els))
            el = lst[t]
            t += l
            yield el
    non_consec_turns = list(ranges(turn_resis))
    non_consec_betas = list(ranges(beta_resis))
    classification_status = []
    #print(turn_resis, non_consec_turns)
    #print(beta_resis, non_consec_betas)
    if len(non_consec_turns) >= 2: # if there are 2 or more non-consecutive turns:
        for t in range(len(non_consec_turns)-1): # for each turn,
            for i in range(len(non_consec_betas)): # for each non-consecutive beta sheet
                #print(f'       {non_consec_turns[1+(1*t)] - non_consec_turns[(1*t)]}')
                if non_consec_turns[1+(1*t)] - non_consec_turns[(1*t)] < 9 and non_consec_turns[1+(1*t)] - non_consec_turns[(1*t)] > 2: # if the start of the 1st turn is at most 9 resis away from the start of the 2nd turn:
                    #print(f'\t    {non_consec_turns[(1*t)]} |     | {non_consec_turns[1+(1*t)]}')
                    if non_consec_betas[i] >= non_consec_turns[(1*t)] and non_consec_betas[i] <= non_consec_turns[1+(1*t)]: # 
                        #print(f'\t    {non_consec_turns[(1*t)]} | {non_consec_betas[i]}     | {non_consec_turns[1+(1*t)]}')
                        classification_status.append("HH")

    if len(classification_status) > 0:
        return 1
    if len(classification_status) == 0:
        return 0

# check if it is a hammerhead or not according to the turn-beta-turn pattern:
def combine_all_dist(betas, turns, hbonds, distances,beta_sheet_distance_cutoff):
    h3_ind = []
    h4_ind = []
    for i in range(len(hbonds)):
        if hbonds[i] == 4:
            h4_ind.append(i)
            h4_ind.append(i+1)
            h4_ind.append(i+2)
            h4_ind.append(i+3)
            h4_ind.append(i+4)
        if hbonds[i] == 3:
            h3_ind.append(i)
            h3_ind.append(i+1)
            h3_ind.append(i+2)
            h3_ind.append(i+3)
    hbonds_arr = np.array(hbonds)
    hbonds_arr[(h4_ind)] = 4
    hbonds_arr[(h3_ind)] = 3
    h = hbonds_arr
    t = []
    b = []
    for i in range(len(h)):
        if i in turns or h[i] == 3 or h[i] == 4:
            t.append('T')
        else:
            t.append('N')
    for i in range(len(hbonds)):
        if i in betas:
            if distances[i] >= beta_sheet_distance_cutoff:
                b.append('D') # beta AND away from HC
            if distances[i] < beta_sheet_distance_cutoff:
                b.append('B') # just beta
        if i not in betas:
            b.append('N')
            
    combined_TB = []
    for i in range(len(t)):
        if b[i] in 'B':
            if t[i] in 'T':
                combined_TB.append('Z') # Z = both a turn and a beta
            if t[i] not in 'T':
                combined_TB.append('B') # B = just a beta sheet, not far from HC 
        if b[i] in 'D':
            if t[i] in 'T':
                combined_TB.append('F')
            if t[i] not in 'T':
                combined_TB.append('D')
        if b[i] not in 'B' and b[i] not in 'D' and t[i] in 'T':
            combined_TB.append('T')
        if b[i] not in 'B' and b[i] not in 'D' and t[i] not in 'T':
            combined_TB.append('N') # N = neither a turn not a beta            
    return combined_TB


def plot_turn_beta_distance(csv_file, name, combined_TB, beta_sheet_distance_cutoff):
    classification_of_hh = check_for_pattern(combined_TB, csv_file['distances'], beta_sheet_distance_cutoff)
    scatter_beta = []
    scatter_turn = []
    scatter_neither = []
    scatter_all = []
    for i in range(len(combined_TB)):
        if combined_TB[i] in 'BDZF':
            scatter_beta.append(1) # beta
        else:
            scatter_beta.append(-1) 
        if combined_TB[i] in 'TZF':
            scatter_turn.append(2) # turn
        else:
            scatter_turn.append(-1)
        if combined_TB[i] in 'N': 
            scatter_neither.append(0) # neither
        else:
            scatter_neither.append(-1)

    beta_colors = []
    for i in range(len(scatter_beta)):
        if scatter_beta[i] == 1 and csv_file['distances'][i] >= beta_sheet_distance_cutoff:
            beta_colors.append('orange')
        if scatter_beta[i] == 1 and csv_file['distances'][i] < beta_sheet_distance_cutoff:
            beta_colors.append('yellow')
        if scatter_beta[i] != 1:
            beta_colors.append('k')

    plt.figure(figsize=[9,2.5])
    if classification_of_hh == 1:
        plt.title(f"{name}\norange: beta + distance >= {beta_sheet_distance_cutoff} A | including type III turn", fontsize=14)
    plt.scatter(np.arange(len(csv_file['ResidueType'])), scatter_turn,    c='blue',   edgecolor='k',    s=200, alpha=1)
    plt.scatter(np.arange(len(csv_file['ResidueType'])), scatter_beta,    c=beta_colors, edgecolor='k', s=200, alpha=1)
    plt.scatter(np.arange(len(csv_file['ResidueType'])), scatter_neither, c='gray',   edgecolor='k',    s=200, alpha=1)
    plt.yticks([0,1,2], ['Neither', 'Beta', 'Turn'], fontsize=16)
    plt.xticks(np.arange(len(csv_file['ResidueType'])), csv_file['ResidueType'], fontsize=16)
    plt.ylim(-0.5,2.5)
    plt.xlim(-0.5,len(csv_file['ResidueType'])-0.5)
    plt.tight_layout()
    plt.savefig(f'plots/{name}.png')
    plt.close()

def classify_turn_beta_distance(csv_file, name, combined_TB, beta_sheet_distance_cutoff):
    classification_of_hh = check_for_pattern(combined_TB, csv_file['distances'], beta_sheet_distance_cutoff)
    return classification_of_hh
