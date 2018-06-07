## PURPOSE: get reads for certain motifs across certain tumors
## INPUT: manifest data 							   all-tumor-manifest.csv
## 		  collapsed fastq files 	sample.converted.unpaired.fastq.collapsed
## OUTPUT: table containing reads for specific motif across samples 	motif.tumor.common.reads.fastq.collapsed.summary.tsv
import os
import os.path
import numpy as np
import pandas as pd
import collections
import subprocess
from pathlib import Path
import time

base_path = '/media/user/2TB (MAC)/Susanna/'
collapsed_ext = '.converted.unpaired.fastq.collapsed'

manifest_file = base_path + 'all-tumor-manifest.csv'
manifest_data =pd.read_csv(manifest_file, header='infer', sep=',')
'''
file = base_path + 'TARGET/TARGET-manifest.csv'
data = pd.read_csv(file, header='infer', sep=',')
data['DISEASE.ABBV'] = 'TARGET'

manifest_data = pd.concat([manifest_data, data])
print(manifest_data.shape)
manifest_data.to_csv(manifest_file, sep=',', index=False)
'''

def getCollapsedFastqDataframe(file):
    df = pd.read_table(file, header=None, delim_whitespace=True)
    df = df.dropna(axis=1, how='all')
    sample = file.split('/')
    sample = sample[len(sample)-1]
    sample = sample.split('.')[0]
    df.columns = ['READS', 'SEQUENCE']
    return df 

def getManifestID(name, tumor):
	id = manifest_data.loc[(manifest_data['DISEASE.ABBV']==tumor) & (manifest_data['NAME']==name)]['ID']
	id = id.tolist()[0]
	id = str(id)
	return str(id)


motifs = ['TGGTTATCTAGCT', 'TTATCAGACTGAT']
mirnas = ['hsa-miR-9-5p-1-2-3', 'hsa-miR-21-5p']

i = 1
motif = motifs[i]
mirna = mirnas[i]


for subdir, dirs, files in os.walk(base_path):
	if '/collapsed_fastq' in subdir:

		folders = subdir.split('/')
		tumor = folders[len(folders)-2]
		if tumor in ['THCA', 'STAD', 'SKCM', 'PCPG']:
			continue
		print(tumor)
		
		summary_file = base_path + 'motif_reads/' + mirna + '/' + motif  +  '.' + tumor + '.common.reads.fastq.collapsed.summary.tsv'
		if Path(summary_file).exists():
			summary_data = pd.read_table(summary_file, header='infer', sep='\t')
		else:
			print('SUMMARY FILE NOT FOUND')
			summary_data = pd.DataFrame({'SEQUENCE':[]})
		matched_ids = list(summary_data)
		common_seqs = list(summary_data['SEQUENCE'])

		total_time_start = time.time()

		for f in os.listdir(subdir):

			time_start = time.time()
			if f[0] == '.':
				break

			patient = f.split('.')[0]
			id = getManifestID(patient, tumor)
			if id not in matched_ids:

				matched_ids.append(id)
				
				if Path(summary_file).exists():
					summary_data = pd.read_table(summary_file, header='infer', sep='\t')
				else:
					print('SUMMARY FILE NOT FOUND')
					summary_data = pd.DataFrame({'SEQUENCE':[]})
				matched_ids = list(summary_data)
				common_seqs = list(summary_data['SEQUENCE'])
				#matched_seq = list(summary_data['SEQUENCE'])
				summary_data = None

				collapsed_file = subdir+'/'+f		
				collapsed_data = getCollapsedFastqDataframe(collapsed_file)
				#print(collapsed_data.shape[0])
				if len(common_seqs) > 0:
					collapsed_data = collapsed_data[collapsed_data.SEQUENCE.isin(common_seqs)]
				num_rows = collapsed_data.shape[0]
				#print(collapsed_data.shape[0])
				collapsed_data.columns = [str(id), 'SEQUENCE']
				match_collapsed_data = collapsed_data #pd.DataFrame(columns = ['READS', 'SEQUENCE'])

				match_collapsed_data.columns = [str(id), 'SEQUENCE']

				if Path(summary_file).exists():
					summary_data = pd.read_table(summary_file, header='infer', sep='\t')
					summary_data = pd.merge(summary_data, match_collapsed_data, on='SEQUENCE', sort=False, how='inner')
				else:
					summary_data = match_collapsed_data


				summary_data.to_csv(summary_file, sep='\t', index=False) 
				summary_data = pd.DataFrame({'SEQUENCE':[]})

				time_end = time.time()

				#print('TUMOR: ' + tumor + ' SAMPLE: ' + str(patient) + ' TOTAL TIME: ' + str((time_end-time_start)/60) + ' ROWS: ' + str(num_rows))
		
		summary_data = pd.read_table(summary_file, header='infer', sep='\t')
		match_summary_data = summary_data.copy()
		for index, row in summary_data.iterrows():
			sequence = str(row['SEQUENCE'])
			if motif not in sequence:
				match_summary_data = match_summary_data[match_summary_data.SEQUENCE != sequence]
		match_summary_data.to_csv(summary_file, sep='\t', index=False) 
		total_time_end = time.time()
		print('TOTAl TUMOR TIME: ' + str(total_time_end-total_time_start)) 
