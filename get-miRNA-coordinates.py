
## PURPOSE: reformat miRNA coordinates from miRBase21 into usable table
## INPUT: miRBase summary data                   miRBase21-master.tsv
##        miRNA coordinates from miRBase        miRNA-coordinates.txt
## OUTPUT: table of miRNA coordinates     miRNA-coordinates-final.csv

import re
import csv
import requests
from lxml import html
import xml.etree.ElementTree
import urllib2
import xmltodict

coordinate_path = '/Users/chens22/Documents/miRNA Stuff/miRNA-coordinates.txt'
summary_path = '/Users/chens22/Documents/miRNA Stuff/RNA Project/miRBase21-master.tsv'
extended_path = '/Users/chens22/Documents/miRNA Stuff/miRNA-coordinates-final.csv'

chromosome = []
x = []
y = []
strand = []
accession = []
mirna = []
sequence = []

def get_reverse_complement(sequence):
        for i in range(0, len(sequence)):
                base = sequence[i]
                if base == 'A':
                        sequence = list(sequence)
                        sequence[i] = 'T'
                        sequence = ''.join(sequence)
                elif base == 'T':
                        sequence = list(sequence)
                        sequence[i] = 'A'
                        sequence = ''.join(sequence)
                elif base == 'C':
                        sequence = list(sequence)
                        sequence[i] = 'G'
                        sequence = ''.join(sequence)                        
                elif base == 'G':
                        sequence = list(sequence)
                        sequence[i] = 'C'
                        sequence = ''.join(sequence)   
        sequence = sequence[::-1]
        return sequence


def get_sequence(chromosome, x, y):
        page = requests.get('http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=' + chromosome + ':' + str(x) + ',' + str(y))
        content = page.content
        
        content = xmltodict.parse(content)
        sequence = content['DASDNA']['SEQUENCE']['DNA']['#text']
        sequence = sequence.replace('\n', '')
        sequence = sequence.replace('\t', '')
        sequence = sequence.upper()
        
        return sequence


sequence_dict = {}
with open(summary_path, 'rb') as csvfile:

        reader = csv.reader(csvfile)
        next(reader, None)
        counter = 0

        for row in reader:
                curr_acc = row[0]
                curr_seq = row[2]
                sequence_dict[curr_acc] = curr_seq.replace('U', 'T')


with open(coordinate_path, 'rb') as f:
        
        data = f.readlines()
        
        for line in data[4:]:
                temp_line = re.split(';|,|\t', line)
                if 'miRNA_primary_transcript' in temp_line:
                        curr_x = int(temp_line[3])-100
                        curr_y = int(temp_line[4])+100
                        curr_chromosome = temp_line[0]
                        curr_strand = temp_line[6]
                        curr_acc = temp_line[8][3:]
                        curr_mirna = temp_line[10][5:len(temp_line[10])-2]
                        
                        chromosome.append(curr_chromosome)
                        x.append(curr_x)
                        y.append(curr_y)
                        strand.append(curr_strand)
                        accession.append(curr_acc)
                        mirna.append(curr_mirna)
                        
                        curr_sequence = get_sequence(curr_chromosome, curr_x, curr_y)
                        if curr_strand in '-':
                                curr_sequence = get_reverse_complement(curr_sequence)
                        
                        old_seq = sequence_dict[curr_acc]
                        if old_seq not in curr_sequence:
                                print 'NOT FOUND: ' + curr_acc

                        sequence.append(curr_sequence)

                
'''
with open(extended_path, 'rb') as csvfile:

        reader = csv.reader(csvfile)
        next(reader, None)
        for row in reader:
                seq = row[6]
                seq = seq.replace('U', 'T')
                extended_seq = row[7]
                print seq
                if seq in extended_seq:
                        print 'FOUND'
                else:
                        print 'NOT FOUND:' + row[4]


'''

with open('/Users/chens22/Documents/miRNA Stuff/miRNA-coordinates-final.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
        headers = ['CHROMOSOME', 'X.COORDINATE', 'Y.COORDINATE', 'STRAND', 'ACCESSION', 'MIRNA', 'SEQUENCE']
        writer.writerow(headers)
        for i in range(0, len(chromosome)):
                row = [chromosome[i], x[i], y[i], strand[i], accession[i], mirna[i], sequence[i]]
                writer.writerow(row) 