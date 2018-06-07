
## PURPOSE: generate unique motifs from the middle bases
## INPUT: miRBase data                            miRBase21-master.tsv
## OUTPUT: table containing unique motifs   mirBase21-motifs(1-10).tsv

import csv
import random

mirna_file = '/Users/chens22/Documents/miRNA/miRBase21-results.tsv'

og_mirna_posi = 0
og_pri_seq_posi = 3
og_mature_seq_posi = 4
og_dupli_posi = 9
og_motif_posi = 6
og_n_motif_posi = 7

pri_seq_posi = 0
mature_seq_posi = 1
dupli_posi = 2
motif_posi = 3
n_motif_posi = 4
n_mseq_posi = 5

def getSummaryData(summary_path):
    
    data = {}
    
    with open(summary_path, 'rU') as tsvin:
        
        tsvin = csv.reader(tsvin, delimiter='\t')
        next(tsvin, None)
        counter = 0
        
        for row in tsvin:

            curr_mirna = row[og_mirna_posi]
            curr_pri_seq = row[og_pri_seq_posi].replace('U', 'T')
            curr_mature_seq = row[og_mature_seq_posi]
            curr_dupli = row[og_dupli_posi]
            #curr_motif = ''#row[og_motif_posi]
            curr_n_motif = row[og_n_motif_posi]
            
            data[curr_mirna] = [curr_pri_seq, curr_mature_seq, curr_dupli, [], [curr_n_motif], '']

    return data


def getSequences(data):
    
    sequences=[]
    
    for m in data.keys():
        sequences+=data[m][pri_seq_posi].split(',')
        
    return sequences


def getMotif(seq, motif_len): #favor first half
    
    return seq[:motif_len]


def isDupli(motif, sequences):
    for seq in sequences:
        if motif in seq:
            return True
    return False


def getMotifs(mature_seq, sequences):
    motif_len = 13

    motifs = []
    start_motif = 0
    
    while start_motif < 1:
        
        curr_motif = getMotif(mature_seq[start_motif:], 13)

        if not isDupli(curr_motif, sequences):
            motifs.append(curr_motif)
            
        start_motif+=1
            
    return motifs


def removeDuplimotifs(motif_dict):
    
    for mirna in motif_dict.keys():
        
        curr_dupli = motif_dict[mirna][dupli_posi]
        
        if curr_dupli in 'TRUE':
            motif_dict[mirna][motif_posi] = ''
            
    return motif_dict


def isNDuplimotif(motif, sequences):
    motif_len = len(motif)
    
    for seq in sequences:
        
        for i in range(0, len(seq)-motif_len+1):
            temp_seq = seq[i:i+motif_len]
            
            for m in range(0, motif_len):

                if motif[m] in 'N':
                    temp_seq = temp_seq[:m] + 'N' + temp_seq[m+1:]
            
            if temp_seq in motif or motif in temp_seq:
                return True
    
    return False


def getNMotif(motif, sequences):

    motif_len = len(motif)
    
    temp_motif = motif[:]
    result_motif = ''
    
    for x in range(1):
        positions = random.sample(range(0, motif_len), motif_len)
        for i in positions:
            temp_motif = temp_motif[:i] + 'N' + temp_motif[i+1:]
                
            if not isNDuplimotif(temp_motif, sequences):
                result_motif = temp_motif[:]
            else:
                temp_motif = motif[:]

    return result_motif


def getNMatureSeq(mseq, n_motifs):
    
    n_mseq = mseq[:]
    
    for motif in n_motifs:
        
        motif_len = len(motif)

        for i in range(0, len(mseq)-motif_len+1):
            temp_mseq = mseq[i:i+motif_len]
            
            for m in range(0, motif_len):
    
                if motif[m] in 'N':
                    temp_mseq = temp_mseq[:m] + 'N' + temp_mseq[m+1:]
            
            if temp_mseq in motif or motif in temp_mseq:
                temp_mseq = mseq[:i] + motif + mseq[i+motif_len:]
                break
        
        for x in range(0, len(n_mseq)):
            if temp_mseq[x] in 'N':
                n_mseq = n_mseq[:x] + 'N' + n_mseq[x+1:]
                
    return n_mseq


def getFastaFile(data):
    
    with open('/Users/chens22/Documents/miRNA/mirBase21-motifs-13(1-10).fa', 'w') as f:
        count = 0
        for mirna in data.keys():
            
            motifs = data[mirna][motif_posi]
            consensus = data[mirna][mature_seq_posi]
            dupli = data[mirna][dupli_posi]
            
            #for m in motifs:
            
            if len(motifs)>0 and dupli in 'FALSE':
                count+=1
                f.write('>' + mirna + '{0: >20}'.format(motifs[0]) + '\n' + consensus + '\n')
        print(count)
        

data = getSummaryData(mirna_file) 
pri_sequences = getSequences(data)

for mirna in data.keys():
    print(mirna)
    curr_mature_seq = data[mirna][mature_seq_posi]
    curr_pri_seq = data[mirna][pri_seq_posi].split(',')
    
    temp_pri_sequences = list(set([x for x in pri_sequences if x not in curr_pri_seq]))
    
    data[mirna][motif_posi] += getMotifs(curr_mature_seq, temp_pri_sequences)
    data[mirna][motif_posi] = list(set(data[mirna][motif_posi]))
    motifs = data[mirna][motif_posi]
data = removeDuplimotifs(data)
getFastaFile(data)
    #for m in motifs:
    #    n_motif = getNMotif(m, temp_pri_sequences)
    #    print(n_motif)
    #    data[mirna][n_motif_posi].append(n_motif)
        
    #data[mirna][n_motif_posi] = list(set(data[mirna][n_motif_posi]))
    #data[mirna][n_mseq_posi] = getNMatureSeq(curr_mature_seq, data[mirna][n_motif_posi])
    
    #print(getNMatureSeq(curr_mature_seq, data[mirna][n_motif_posi]))
    #print(curr_mature_seq)
    



with open('/Users/chens22/Documents/miRNA/mirBase21-motifs(1-10).tsv', 'w') as f:
    headers = ['MIRNA', 'SEQUENCE', 'CONSENSUS', 'DUPLIMOTIF', 'MOTIF.13']#, 'N.MOTIF', 'N.CONSENSUS'] 
    f.write('\t'.join(headers) + '\n')
    for mirna in data.keys():
        row = [mirna]
        
        for i in range(0, len(headers)-1):
            curr_element = data[mirna][i]
            row.append(str(curr_element))

        f.write('\t'.join(row) + '\n') 

'''

def get_in_posi(n):

    posi = []
    
    for i in range(0,n/2):
        posi.append(i)
        posi.append(n-i-1)
        
    if n%2==1:
        posi.append(n/2)
    
    return posi

'''





