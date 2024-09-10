# 06102024
# converts a numpy array of full HC sequences to individual fasta sequences, each named sequence_n.fasta


#LC seq: https://www.ncbi.nlm.nih.gov/nuccore/GU272046.1

LC_seq = 'QSALTQPASVSGSPGQSITISCNGTSNDVGGYESVSWYQQHPGKAPKVVIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEGDYYCKSLTSTRRRVFGTGTKLTVL'
import numpy as np
full_seqs = np.load('PG9_noVgene_cdrh3_anylen_motif_4__full_seqs.npy')

for i in range(len(full_seqs)):
    f = open(f"a3m_files/sequence_{i+1}.a3m", "w")
    f.write(f"#{len(full_seqs[i])},108	1,1\n")
    
    f.write(">101	102\n")
    f.write(f"{full_seqs[i]}{LC_seq}\n")
    
    f.write(">101\n")
    LC_dashes = "-"*(len(LC_seq))
    f.write(f"{full_seqs[i]}{LC_dashes}\n")
    
    f.write(">102\n")
    HC_dashes = "-"*(len(full_seqs[i]))
    f.write(f"{HC_dashes}{LC_seq}\n")
    f.close()
