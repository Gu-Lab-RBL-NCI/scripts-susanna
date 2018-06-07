## PURPOSE: get 2d struct for each added base to miRNA sequence to determine changes in structure
## INPUT: miRBase data with extended sequences                                      miRBase21-master.tsv
## OUTPUT: table containing 2d struct with each lengthen in sequence    primirna-energy-extended-seq.tsv

import os
import csv
from subprocess import call
import tempfile

mirna_file = '/Users/chens22/Documents/miRNA/miRBase21-master.tsv'

primirna = 'hsa-mir-9-1'

def getExtendedSequence(miRBase_file, mirna, coordinates):
    with open(mirna_file,'rU') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            curr_primirna = row[2]
            if(curr_primirna in mirna):
                pri_seq = row[3].replace('T', 'U') 
                extended_seq = row[20].replace('T', 'U')
                start = extended_seq.find(pri_seq)
                end = start + len(pri_seq)
                start = start + coordinates[0]
                end = end + coordinates[1]
                return(extended_seq[start:end])
            

def getSecondaryStructure(sequence, primirna):
    tempfile = '/Users/chens22/Documents/miRNA/scripts/temp.fa'
    with open(tempfile, 'w') as text_file:
        
        text_file.write('>' + primirna + '\n')
        text_file.write(sequence + '\n')
        text_file.close()
        
    structfile = 'temp.txt'
    open(structfile, 'w').close()
    os.system('RNAfold --infile=\'' + tempfile + '\' --outfile=\'' + structfile + '\'')
    
    with open(structfile, 'r') as f:
        
        for line in f:
            if('>' in line):
                mirna = line.replace('>', '').replace('\n', '')
                next(f)
                data = next(f)
                data = data.split(' ')
                struct = data[0]
                energy = data[1].replace('(', '').replace(')', '').replace('\n', '')
                print(data)
                return([struct, energy]) 
        
#print(getSecondaryStructure(getPrimarySequence(mirna_file, primirna), primirna))
filename = '/Users/chens22/Documents/miRNA/' + primirna + '-energy-extended-seq.tsv'
with open(filename, 'w') as tsvout:
    
    tsvout.write('MIRNA\tSEQUENCE\tCOORDINATES\tSECONDARY.STRUCTURE\tTOTAL.ENERGY\tENERGY.PER.BASE')
    
    coordinates = (0,0)
    for i in range(0, 20):
        coordinates = (-i-1, i)
        extended_seq = getExtendedSequence(mirna_file, primirna, coordinates)
        secondary_structure = getSecondaryStructure(extended_seq, primirna)
        struct = secondary_structure[0]
        energy = secondary_structure[1]
        energy_per_base = float(energy)/len(extended_seq)
        tsvout.write('\n' + primirna + '\t' + extended_seq + '\t' + str(coordinates) + '\t' + struct + '\t' + energy + '\t' + str(energy_per_base))

        coordinates = (-i-1, i+1)
        extended_seq = getExtendedSequence(mirna_file, primirna, coordinates)
        secondary_structure = getSecondaryStructure(extended_seq, primirna)
        struct = secondary_structure[0]
        energy = secondary_structure[1]
        energy_per_base = float(energy)/len(extended_seq)
        tsvout.write('\n' + primirna + '\t' + extended_seq + '\t' + str(coordinates) + '\t' + struct + '\t' + energy + '\t' + str(energy_per_base))
    
        