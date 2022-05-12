from Bio.PDB.PDBParser import PDBParser
import os 
import glob
import re 

drug_name = 'Ibuprofen'
CWD = os.getcwd() + f'/data/{drug_name}/'
pattern = '/([aA0-zZ9]*).pdb'
pdb_records = glob.glob(CWD+'*.pdb')

for pdb_record in pdb_records:

    ID = re.findall(pattern, pdb_record)[0]

    # You can use a dict to convert three letter code to one letter code
    t3t1 = {'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    }

    # run parser
    parser = PDBParser(QUIET=True)
    structure_2 = parser.get_structure(f'{ID}',pdb_record)    

    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    for model in structure_2:
        for chain in model:
            seq = []
            residue_n = 0
            residue_ibu = 0
            for residue in chain:
                if residue.resname in t3t1:
                    seq.append(t3t1[residue.resname])
                    #print('>OXIDOREDUCTASE\n',''.join(seq))
                
                #else: 
                    #print(f'Model: {model}\nChain: {chain}\nResidue_n: {residue_n}\nResidue: {residue.resname}')
                    
                
                residue_n += 1

    protein_fasta = ''.join(seq)




