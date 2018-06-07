## PURPOSE: SUMMARIZE ALL READS FOR A TUMOR
## INPUT: collapsed fastq files      sample.converted.unpaired.fastq.collapsed
## OUTPUT: summary collapsed fastq file     



from os import listdir
from os.path import isfile, join
import numpy
import pandas as pd

base_path = '/Users/chens22/Documents/miRNA/'


tumors = ['ACC']
for t in tumors:
    fastq_path = base_path + t + '/collapsed_fastq/'
    collapsed_files = [f for f in listdir(fastq_path) if (isfile(join(fastq_path, f)) and not f.startswith('.'))]    
    
    for c in collapsed_files:
        if '.1.' not in c:
            df = pd.read_table(base_path + t + '/collapsed_fastq/' + c, header=None, delim_whitespace=True)
            df = df.dropna(axis=1, how='all')
            
        summary = {}
        for index, row in df.iterrows():
            curr_len = len(summary)
            if curr_len%10000==0:
                print curr_len
            seq = row[1]
            reads = row[0]
            if seq not in (summary.keys()):
                summary[seq] = reads
            else:
                summary[seq] = summary[seq] + reads
        
        sample = c.split('.')[0]
        filename = base_path + t + '/collapsed_fastq/' + sample + '.collapsed.fastq.summary.tsv'
        
        with open(filename, 'w') as f:
            f.write('SEQUENCE\tREADS')
            for seq in list(summary.keys()):
                f.write(seq + '\t', + summary[seq])
            
            