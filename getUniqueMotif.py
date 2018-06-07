
## PURPOSE: get unique motif for indistinguishable miRNA using IUPAC code with shifting distance starting from center
## INPUT: miRBase data           miRBase21-master.tsv
## OUTPUT: table containing motifs  unique_motifs.tsv

import numpy as np
import pandas as pd
import difflib

mirna_file = '/Users/chens22/Documents/miRBase21/miRBase21-master.tsv'
mirna_data = pd.read_table(mirna_file, header='infer', sep='\t')

def getMatureSeqDist(motif, sequences):
    dupli_seq = 0
    dist = 0
    for s in sequences:
        if motif in s:
            dupli_seq += 1
        else:
            for d in difflib.ndiff(motif, s):
                if d[0]=='-' or d[0]=='+':
                    dist += 1  
    if len(sequences)-dupli_seq < 1:
        return 0
    return (dist/(len(sequences)-dupli_seq))-dupli_seq

def inMatureSeq(motif, sequences):
    for s in sequences:
        if motif in s:
            return True
    return False

def translateIUPAC(seq):
    codes = {'N':['A', 'G', 'C', 'T'],
             'V':['A', 'G', 'C'],
             'H':['A', 'T', 'C'],
             'D':['A', 'G', 'T'],
             'B':['T', 'G', 'C'],
             'M':['A', 'C'],
             'K':['G', 'T'],
             'W':['A', 'T'],
             'S':['G', 'C'],
             'Y':['T', 'C'],
             'R':['A', 'G'],
             }    
    

def inIUPACMatureSeq(motif, sequences):
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

        
def getMotif(seq, mid, length):
    start = mid - int(length/2)
    end = mid + int(round(length/2)) + 1
    motif = seq[start:end]
    return motif

def getAltNumSeq(n):
    result = []
    for i in range(0, round(n/2)):
        result.append(i)
        result.append(n-i)
    return result

def getDistanceFromCenter(motif, seq):
    start_margin = seq.find(motif)
    end_margin = len(seq) - (start_margin + len(motif))
    return end_margin - start_margin

def getIUPAC(base):
    codes = {'N':['A', 'G', 'C', 'T'],
             'V':['A', 'G', 'C'],
             'H':['A', 'T', 'C'],
             'D':['A', 'G', 'T'],
             'B':['T', 'G', 'C'],
             'M':['A', 'C'],
             'K':['G', 'T'],
             'W':['A', 'T'],
             'S':['G', 'C'],
             'Y':['T', 'C'],
             'R':['A', 'G'],
            }
    result = []
    for key in codes:
        values = codes[key]
        if base in values:
            result.append(key)
    return result

def chop_motif(motif):
    n_pos = str.find(motif, 'N')
    if n_pos==-1:
        return([motif])
    n_num = 0
    i = n_pos
    while i<len(motif) and motif[i]=='N':
        i+=1
        n_num+=1
    return([motif[:n_pos], n_num] + chop_motif(motif[n_pos+n_num:]))

print(chop_motif('NNNTAGAGGGAAGCGCTTTNNN'))

'''
dupli_mirna = list(set(mirna_data.loc[mirna_data['DUPLIMOTIF']==True]['MIRNA']))
final_motifs = {}
for d in dupli_mirna:

    print(d)
    curr_seq = str(list(set(mirna_data.loc[mirna_data['MIRNA']==d]['SEQUENCE']))[0])
    print(curr_seq)
    sequences = list(set(mirna_data.loc[mirna_data['MIRNA']!=d]['SEQUENCE']))

    mid = int(len(curr_seq)/2)
    bound = 3
    max_dist = 0
    motif = getMotif(curr_seq, mid, 13)
    most_unique_motif = motif
    for i in range(-bound, bound):
        dist = getMatureSeqDist(motif, sequences)
        if dist > max_dist:
            max_dist = dist
            most_unique_motif = motif
            shift = i
   
    motif = most_unique_motif
    print(motif)
    
    start_motif = curr_seq.find(motif)
    end_motif = start_motif + len(motif)
    
    while inMatureSeq(motif, sequences):
        if start_motif > 0:
            start_motif-=1
        motif = curr_seq[start_motif:end_motif]
        if end_motif < len(curr_seq):
            if inMatureSeq(motif, sequences):
                end_motif+=1
        else:
            break
        motif = curr_seq[start_motif:end_motif]
    print(motif)
        
    start_motif = curr_seq.find(motif)
    end_motif = start_motif + len(motif)    
        
    if len(motif) > 13:
        for i in range(0, round(len(motif)/2)):
            temp_start_motif = start_motif + 1
            temp_motif = curr_seq[temp_start_motif:end_motif]
            if not inMatureSeq(temp_motif, sequences) and len(temp_motif)>=13:
                start_motif = temp_start_motif
            
            temp_end_motif = end_motif - 1
            temp_motif = curr_seq[start_motif:temp_end_motif]
            if not inMatureSeq(temp_motif, sequences) and len(temp_motif)>=13:
                end_motif = temp_end_motif

        motif = curr_seq[start_motif:end_motif]
            
    print(motif)
    
    dist_cent = getDistanceFromCenter(motif, curr_seq)
    print('DISTANCE FROM CENTER: ' + str(dist_cent))
    jump = 1
    if dist_cent < 0:
        jump = -1
    for i in range(0, dist_cent, jump):
        temp_motif = curr_seq[start_motif + i:end_motif+i]
        if not inMatureSeq(temp_motif, sequences):
            motif = temp_motif

    
    if dist_cent < 0:
        start_motif += dist_cent
    else:
        end_motif += dist_cent
    
    motif = curr_seq[start_motif:end_motif]
    print(motif)
    
    dist_cent = getDistanceFromCenter(motif, curr_seq)
    print('DISTANCE FROM CENTER: ' + str(dist_cent))   
     
    
    if len(motif) > 13:
        for i in range(1, int((len(motif)-13)/2)):
            temp_motif = 'N'*i + motif[i:]
            if inIUPACMatureSeq(temp_motif, sequences):
                motif = temp_motif
            
            temp_motif = motif[:len(motif)-i] + 'N'*i
            if inIUPACMatureSeq(temp_motif, sequences):
                motif = temp_motif 
        
        start = 0
        end = len(motif)
        for i in range(0, round(len(motif)/2)):
            temp_motif = motif[i:end]
            if not inIUPACMatureSeq(temp_motif, sequences) and len(temp_motif)>=13:
                start = i
            temp_motif = motif[start:len(motif)-i]
            if not inIUPACMatureSeq(temp_motif, sequences) and len(temp_motif)>=13:
                end = len(motif)-i
        motif = motif[start:end]
            
    print(motif)
    final_motifs[d] = [motif, curr_seq, dist_cent]
    
with open('/Users/chens22/Documents/mirBase21/unique-duplimotifs.tsv', 'w') as f:
    headers = ['MIRNA', 'MOTIF', 'SEQUENCE', 'DISTANCE.FROM.CENTER']
    f.write('\t'.join(headers) + '\n')
    for m in final_motifs.keys():
        row = [m]
        
        for i in range(0, len(headers)-1):
            curr_element = final_motifs[m][i]
            row.append(str(curr_element))

        f.write('\t'.join(row) + '\n') 
    
'''
        
        

        
        
        
   
    
    