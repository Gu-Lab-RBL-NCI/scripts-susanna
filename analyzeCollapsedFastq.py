
## PURPOSE: return file containing reads from all collapsed fastq files in tumor and/or sequence overlap between two samples
## INPUT: collapsed fastq files     sample.converted.unpaired.fastq.collapsed
## OUTPUT: collapsed fastq summary  tumor.converted.unpaired.fastq.collapsed.summary.tsv

from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd

base_path = '/Users/chens22/Documents/miRNA/'


def getNRowDataframe(df):
    nrow = df.shape[0]
    return nrow

def getSequenceOverlap(seq_list1, seq_list2):
    
    if not isinstance(seq_list1, pd.DataFrame):
        seq_list1 = pd.DataFrame(seq_list1)
    if not isinstance(seq_list2, pd.DataFrame):
        seq_list2 = pd.DataFrame(seq_list2)     
    
    max_len = getNRowDataframe(seq_list1)
    len2 = getNRowDataframe(seq_list2)
    if max_len <= len2:
        max_len = len2     
    
    seq_list = seq_list1.append(seq_list2)
    seq_list = seq_list[seq_list.duplicated(keep=False)]
    seq_list = seq_list.drop_duplicates()
    
    return float(getNRowDataframe(seq_list)/float(max_len))


def getFastqDataframe(file):
    df = pd.read_table(file, header=None, delim_whitespace=True)
    df = df.dropna(axis=1, how='all')
    sample = file.split('/')
    sample = sample[len(sample)-1]
    sample = sample.split('.')[0]
    df.columns = [str(sample) + '.READS', 'SEQUENCE']
    return df 


tumors = ['ACC', 'BRCA', 'CHOL', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
          'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
          'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
#BLCA, COAD, CESC
tumors = ['ACC']


for t in tumors:
    result_file = base_path + t + '/' + t + '.converted.unpaired.fastq.collapsed.summary.tsv'          
    fastq_path = base_path + t + '/collapsed_fastq/'
    collapsed_files = [f for f in listdir(fastq_path) if (isfile(join(fastq_path, f)) and not f.startswith('.'))]    
    
    if len(collapsed_files) > 0:
        print '\n' + t
    
    reads_data = pd.DataFrame({'SEQUENCE': []})
    total_rows = 0
    for c in collapsed_files:
    
        file = base_path + t + '/collapsed_fastq/' + c
        if '.1.' not in c or '.2.' not in c:
            df = getFastqDataframe(file)
            reads_data = pd.merge(reads_data, df, on='SEQUENCE', sort=False, how='outer')
            
            #print "SEQUENCE OVERLAP: " + str(getSequenceOverlap(df1['SEQUENCE'], df2['SEQUENCE']))
        else:
            os.remove(file)
        
        reads_data.to_csv(result_file, sep='\t') 
        reads_data.to_csv(result_file, sep='\t')         
         
