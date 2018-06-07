## PURPOSE: collapse fastq files
## INPUT: fastq files	sample.converted.unpaired.fastq
## OUTPUT: collapsed files	sample.converted.unpaired.fastq.collapsed

import os
import os.path
import collections
import subprocess
'''
def batch_iterator(iterator, batch_size):
	entry = True
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.next()
			except StopIteration:
				entry = None
			if entry is None:
				break
			batch.append(entry)
		if batch:
			yield batch

def getCollapsedFastqDataframe(file):
	df = pd.read_table(file, header=None, delim_whitespace=True)
	df = df.dropna(axis=1, how='all')
	df.columns = ['READS', 'SEQUENCE']
	return df 
'''

base_path = '/Volumes/2TB (MAC)/Susanna/'
collapsed_ext = '.converted.unpaired.fastq.collapsed'

collapsed = {}
for subdir, dirs, files in os.walk(base_path):
	if '/fastq' in subdir:

		folders = subdir.split('/')
		tumor = folders[len(folders)-2]
		print tumor
		if tumor not in ['COAD', 'BLCA']:
			continue
		
		for f in os.listdir(subdir):
			patient = f.split('.')[0]
			fastq_file = subdir+'/'+f		
			outfile = base_path + tumor + '/collapsed_fastq/' + patient + collapsed_ext	
			if not os.path.exists(outfile) or os.path.getsize(outfile) < 1 or 1 == 1:
				print "WORKING ON: " + outfile
				os.system('awk "NR%4==2" "'+ fastq_file + '" | sort -S1900G | uniq -c > "' + outfile + '"')

			