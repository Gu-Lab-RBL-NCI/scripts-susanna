## PURPOSE: download gene and miRNA expression files from CGC
## INPUT: manifest data 	all-tumor-manifest.csv
## OUTPUT: downloaded files
#from __future__ import division
import sevenbridges as sbg
import os
import numpy as np
import pandas as pd
from pathlib import Path
import collections
import sys
from math import ceil
import json
from requests import request
from collections import defaultdict

base_path = '/Users/chens22/Documents/'
collapsed_ext = '.converted.unpaired.fastq.collapsed'
fastq_ext = '.converted.unpaired.fastq'

base_url = 'https://cgc-datasets-api.sbgenomics.com/datasets/target/v0/'

manifest_file = base_path + 'all-tumor-manifest.csv'
manifest_data = pd.read_table(manifest_file, sep=',', header='infer')
#manifest_data = manifest_data.dropna()

type = 'miRNA'
if type in 'miRNA':
	ext = '.mirnas.quantification.txt'#'.FPKM-UQ.txt.gz'#
	folder = '/miRNA_expression/'#'/gene_expression/'#
else:
	ext = '.FPKM-UQ.txt.gz'#
	folder = '/gene_expression/'#


def api_call(path, method='GET', query=None, data=None, token=None):
     
    #base_url = 'https://cgc-api.sbgenomics.com/v2/'

    data = json.dumps(data) if isinstance(data, dict) \
    or isinstance(data,list) else None
               
    headers = {
        'X-SBG-Auth-Token': token,
        'Accept': 'application/json',
        'Content-type': 'application/json',
        'X-RateLimit-Limit': str(1000)
    }
     
    response = request(method, base_url + path, params=query, \
                       data=data, headers=headers)
    response_dict = response.json() if response.json() else {}
 
    if response.status_code / 100 != 2:
        print(response_dict)
       # print('Error Code: %i.' % (response_dict['code']))
        #print(response_dict['more_info'])
        raise Exception('Server responded with status code %s.' \
                        % response.status_code)
    return response_dict



auth_token = 'b5764171da2e446f9fffdf2047eb0035'
os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
os.environ['SB_AUTH_TOKEN'] = auth_token
api = sbg.Api()


tumors = list(set(list(manifest_data['DISEASE.ABBV'])))
tumors = [str(x) for x in tumors]
tumors = [str(x) for x in tumors if x != 'nan']
for tumor in tumors:

	#if tumor in ['TARGET']:
	#	continue
	if tumor not in 'TARGET':
		continue
	diseases = list(set(list(manifest_data.loc[manifest_data['DISEASE.ABBV']==tumor,'DISEASE.TYPE'])))

	diseases = ['']
	for d in diseases:

		disease = d.lower().title().replace('.', ' ').replace('And', 'and').replace('-C', '-c').replace('High ', 'High-')

		disease = 'Rhabdoid Tumor'
		data_type = type + " Expression Quantification"
		if type in 'gene':
			data_type = data_type.title()
		query_body = {
		    "entity": "files",
		    "hasDataType" : data_type,
		    #"hasDataFormat": "TXT",
		    "hasInvestigation": {
		        "hasDiseaseType" : disease
		    }

		}


		if tumor not in 'TARGET':
			disease = ''

		print(tumor)
		total = api_call(method='POST', path ='query/total', token=auth_token, data=query_body)

		if total['total'] == 0:
			print('NO ' + type + ' EXPRESSION FILES FOUND:' + tumor)
			print(disease)
			continue

		print(total)

		files_in_query = []
		 
		loops = int(ceil(total['total']/100))
		#loops = 1
		for ii in range(0,loops):
			files_in_query.append(api_call(method='POST', path ='query/?offset=' + str(100*ii), token=auth_token, data=query_body)) #("query/?offset=100?limit=3")

		file_list = []

		for ii in range(0, loops):
			for f in files_in_query[ii]['_embedded']['files']:
				if ext in f['label'] or f['label'] in ext:

					try:
						curr_file = api.files.get(id = f['id'])
					except:
						if str(os.environ['SB_API_ENDPOINT']) in base_url:
							os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
						else:
							os.environ['SB_API_ENDPOINT'] = base_url
						api = sbg.Api()
						curr_file = api.files.get(id = f['id'])
					file_list.append(curr_file)

		dl_dir = base_path + tumor + '/'
		try:
		    os.stat(dl_dir)
		except:
		    os.mkdir(dl_dir)

		dl_dir = base_path + tumor + folder
		try:
		    os.stat(dl_dir)
		except:
		    os.mkdir(dl_dir)

		dl_dir = base_path + tumor + folder + disease + '/'
		try:
		    os.stat(dl_dir)
		except:
		    os.mkdir(dl_dir)
		 
		if tumor in 'TARGET':
			if 'gene' in folder:
				meta_file = base_path + tumor + '/gene_expression/'  + tumor + '-' + disease + '-' + type + '-expression-manifest.csv'
			else:
				meta_file = base_path + tumor + '/miRNA_expression/' + tumor  + '-' + disease + '-' + type + '-expression-manifest.csv'
		
		else:
			meta_file = base_path + tumor + disease + '/' + tumor + '-' + type + '-expression-manifest.csv'
	
		metadata = defaultdict(list)
		#if not os.path.exists(meta_file):
		for f in file_list:
			try:
				curr_metadata = dict(f.metadata)
			except:
				if str(os.environ['SB_API_ENDPOINT']) in base_url:
					os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
				else:
					os.environ['SB_API_ENDPOINT'] = base_url
				api = sbg.Api()
				try:
					curr_file = api.files.get(id = f.id)
					file_list.remove(f)
					file_list.append(curr_file)
					curr_metadata = dict(curr_file.metadata)
				except:
					if str(os.environ['SB_API_ENDPOINT']) in base_url:
						os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
					else:
						os.environ['SB_API_ENDPOINT'] = base_url
					api = sbg.Api()
					curr_file = api.files.get(id = f.id)
					file_list.remove(f)
					file_list.append(curr_file)
					try:
						curr_metadata = dict(curr_file.metadata)
					except:
						if str(os.environ['SB_API_ENDPOINT']) in base_url:
							os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
						else:
							os.environ['SB_API_ENDPOINT'] = base_url
						api = sbg.Api()
						curr_file = api.files.get(id = f.id)
						curr_metadata = dict(curr_file.metadata)
						'''
						try:
							curr_file = api.files.get(id = f.id)
							file_list.remove(f)
							file_list.append(curr_file)
							curr_metadata = dict(curr_file.metadata)
						except:
							if str(os.environ['SB_API_ENDPOINT']) in base_url:
								os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
							else:
								os.environ['SB_API_ENDPOINT'] = base_url
							api = sbg.Api()
							curr_file = api.files.get(id = f.id)
							file_list.remove(f)
							file_list.append(curr_file)
							try:
								curr_metadata = dict(curr_file.metadata)
							except:
								if str(os.environ['SB_API_ENDPOINT']) in base_url:
									os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
								else:
									os.environ['SB_API_ENDPOINT'] = base_url
								api = sbg.Api()
								try:
									curr_file = api.files.get(id = f.id)
									file_list.remove(f)
									file_list.append(curr_file)
								except:
									if str(os.environ['SB_API_ENDPOINT']) in base_url:
										os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
									else:
										os.environ['SB_API_ENDPOINT'] = base_url
									api = sbg.Api()
									curr_file = api.files.get(id = f.id)
								curr_metadata = dict(curr_file.metadata)
							'''
			
			print(curr_metadata)
			curr_metadata['id'] = f.id
			curr_metadata['name'] = curr_metadata['aliquot_id'] + ext
			if 'gene' not in folder:
				curr_metadata['name'] = curr_metadata['aliquot_id'] + '.mirna.quantification.txt'
			
			for key, value in curr_metadata.items():
				metadata[key].append(value)
		metadata = pd.DataFrame.from_dict(metadata, orient='columns')
		metadata.to_csv(meta_file, sep=',')
'''
		for f in file_list:
			if type in 'gene':
				name = f.name
				curr_file = f
			else:
				curr_file = f
				try:
					name = curr_file.metadata['aliquot_id']
				except:
					os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'

					api = sbg.Api()
					curr_file = api.files.get(id = f.id)
					print(curr_file)
					name = curr_file.metadata['aliquot_id']
				ext = '.mirna.quantification.txt'
			if not os.path.exists(dl_dir + name + ext):
				try:
					curr_file.download(path = (dl_dir + '/' + name + ext))
				except:
					if str(os.environ['SB_API_ENDPOINT']) in base_url:
						os.environ['SB_API_ENDPOINT'] = 'https://cgc-api.sbgenomics.com/v2/'
					else:
						os.environ['SB_API_ENDPOINT'] = base_url
					api = sbg.Api()
					file = api.files.get(id = f.id)
					print((dl_dir + '/' + name + ext))
					file.download(path = (dl_dir + '/' + name + ext))


		#get metadata
		
		
'''
		