## PURPOSE: crop 3d miRNA structures
## INPUT: miRBase data      miRBase21-master.tsv
##        miRNA 3d structures   mirna.pdb
## OUTPUT: crops 3d structures to bases 1-33

from Bio.PDB import PDBParser
from Bio.PDB.PDBParser import PDBParser
from Bio import PDB
from Bio.PDB import PDBIO
import sys
import csv
import Bio.PDB

mirna_file = '/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'

low_fid_mirna = ['hsa-miR-758-5p', 'hsa-miR-153-5p','hsa-miR-383-5p','hsa-miR-382-5p','hsa-miR-183-5p',
                 'hsa-miR-425-5p','hsa-miR-136-5p','hsa-miR-628-5p','hsa-miR-769-5p','hsa-miR-30e-5p',
                 'hsa-miR-98-5p','hsa-miR-138-5p-1-2','hsa-miR-144-5p','hsa-let-7i-5p', 'hsa-miR-191-5p',
                 'hsa-miR-532-5p']
high_fid_mirna = ['hsa-miR-100-5p','hsa-miR-106b-5p','hsa-miR-129-5p-1-2','hsa-miR-589-5p','hsa-miR-125b-5p-1-2',
                  'hsa-miR-30a-5p','hsa-miR-125a-5p','hsa-miR-134-5p','hsa-miR-99b-5p','hsa-miR-330-5p',
                  'hsa-miR-26a-5p-1-2','hsa-miR-21-5p','hsa-let-7b-5p','hsa-miR-1307-5p','hsa-miR-185-5p',
                  'hsa-miR-361-5p','hsa-miR-379-5p','hsa-let-7f-5p-1-2','hsa-let-7e-5p','hsa-miR-340-5p',
                  'hsa-miR-149-5p','hsa-miR-204-5p','hsa-miR-26b-5p','hsa-miR-30d-5p']
def getPrimirna(mirna):
    primirna = []
    with open(mirna_file, 'rU') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')     
        next(tsvin, 0)
        for row in tsvin:
            curr_mirna = row[0]
            curr_primirna =  row[2]
            if curr_mirna in mirna:
                primirna.append(curr_primirna)
    return(primirna)

def getResidueCount(struct):
    res_no = 0
    for model in struct:
        for chain in model:
            for r in chain.get_residues():
                #print('count' + str(r.id[1]))
                if r.id[0] == ' ':
                    res_no+=1
    return(res_no)

def getAlignedStructure(file1, file2, start_id, end_id):
    
    atoms_to_be_aligned = range(start_id, end_id + 1)
    
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    
    # Get the structures
    ref_structure = pdb_parser.get_structure("reference", file1)
    sample_structure = pdb_parser.get_structure("sample", file2)
    
    ref_model    = ref_structure[0]
    sample_model = sample_structure[0]
    
    ref_atoms = []
    sample_atoms = []
    
    for ref_chain in ref_model:
        for ref_res in ref_chain:
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                ref_atoms.append(ref_res)
    
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res)
    
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    
    # Print RMSD:
    print super_imposer.rms
    
    return sample_structure

    
def cropStructure(structure, start_id, end_id):

    res_count = getResidueCount(structure)
    
    resi_ids = []
    for model in structure:
        for chain in model:
            for r in chain.get_residues():
                if r.id[1]>end_id or r.id[1] < start_id: #determines which residues to keep
                    resi_ids.append(r.id)
    
    for model in structure:
        for chain in model:
            for r in resi_ids:
                chain.detach_child(r)

    return structure 








pdb1 = '/Users/chens22/Dropbox/miRNA Analysis/tertiary structure/hsa0mir0758.pdb'
pdb2 = '/Users/chens22/Dropbox/miRNA Analysis/tertiary structure/hsa0mir0758corr.pdb'

struct1 = PDBParser().get_structure('struct1', pdb1)
struct2 = PDBParser().get_structure('struct2', pdb2)

struct1 = cropStructure(struct1, 1, 33)
struct2 = cropStructure(struct2, 1, 33)
parser = PDBParser()
w = PDBIO()
w.set_structure(struct1)
w.save(pdb1) 
w.set_structure(struct2)
w.save(pdb2) 

aligned_struct = getAlignedStructure(pdb1, pdb2, 1, 33)



'''
for m in high_fid_mirna:
    primirna = getPrimirna(m)
    for p in primirna:
        p = p.replace('-', '0')
        pdb = '/Users/chens22/Documents/miRNA/structure/high fidelity tertiary structure/' + p + '.pdb'
        structure = cropStructure(pdb)
        parser = PDBParser()
        w = PDBIO()
        w.set_structure(structure)
        w.save(pdb)           
'''   
