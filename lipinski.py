#!/usr/bin/env python
'''
Chemical data about a molecule.

Molecules are defined by SMILES strings. Can work out logP values, Lipinski's 
rules, etc...

Uses rdkit

SOURCE: https://gist.github.com/fdc4db6d450b66345f46.git
'''

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
import sys
import pandas as pd
import argparse

class SmilesError(Exception): pass

def log_partition_coefficient(smiles):
    '''
    Returns the octanol-water partition coefficient given a molecule SMILES 
    string
    '''
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        raise SmilesError('%s returns a None molecule' % smiles)
        
    return Crippen.MolLogP(mol)
    
def lipinski_trial(smiles):
    '''
    Returns which of Lipinski's rules a molecule has failed, or an empty list
    
    Lipinski's rules are:
    Hydrogen bond donors <= 5
    Hydrogen bond acceptors <= 10
    Molecular weight < 500 daltons
    logP < 5
    '''
    passed = []
    failed = []
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise Exception('%s is not a valid SMILES string' % smiles)
    
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)
    
    failed = []
    
    if num_hdonors > 5:
        failed.append('Over 5 H-bond donors, found %s' % num_hdonors)
    else:
        passed.append('Found %s H-bond donors' % num_hdonors)
        
    if num_hacceptors > 10:
        failed.append('Over 10 H-bond acceptors, found %s' \
        % num_hacceptors)
    else:
        passed.append('Found %s H-bond acceptors' % num_hacceptors)
        
    if mol_weight >= 500:
        failed.append('Molecular weight over 500, calculated %s'\
        % mol_weight)
    else:
        passed.append('Molecular weight: %s' % mol_weight)
        
    if mol_logp >= 5:
        failed.append('Log partition coefficient over 5, calculated %s' \
        % mol_logp)
    else:
        passed.append('Log partition coefficient: %s' % mol_logp)
    
    return passed, failed
    
def lipinski_pass(smiles):
    '''
    Wraps around lipinski trial, but returns a simple pass/fail True/False
    '''
    passed, failed = lipinski_trial(smiles)
    if failed:
        return False
    else:
        return True
'''    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python lipinski.py <input_file.smi>")
        sys.exit(1)

    input_file = sys.argv[1]

    try:
        df = pd.read_csv(input_file)
        if 'smiles' not in df.columns:
            print(f"Warning: The input file does not contain a 'smiles' column. Using the first column as 'smiles'.")
            df.columns = ['smiles'] + df.columns.tolist()[1:]
        smiles_list = df['smiles'].tolist()
    except IOError:
        print(f"Error: Could not read file {input_file}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: The input file is empty")
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Error: The input file is not a valid CSV")
        sys.exit(1)

    for smiles in smiles_list:
        smiles = smiles.strip()
        if smiles:
            result = lipinski_pass(smiles)
            print(f"SMILES: {smiles} - Pass: {result}")
'''

def main():
    parser = argparse.ArgumentParser(description='Process some SMILES strings.')
    parser.add_argument('input_file', type=str, help='Input CSV file containing SMILES strings')
    parser.add_argument('-o', '--output', type=str, default='result_lipinski.csv', help='Output CSV file (default: result_lipinski.csv)')
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file)
        if 'smiles' in df:
            smiles_list = df['smiles'].tolist()
        elif 'Smiles' in df:
            smiles_list = df['Smiles'].tolist()
        elif 'SMILES' in df:
            smiles_list = df['SMILES'].tolist()
        else: 
            print(f"Warning: The input file does not contain a 'smiles' column. Using the first column as 'smiles'.")
            smiles_list = pd.read_csv(args.input_file, header=None).iloc[:, 0].tolist()
    except FileNotFoundError:
        print(f"Error: The file {args.input_file} does not exist")
        sys.exit(1)
    except KeyError:
        print(f"Warning: The input file does not contain a 'smiles' column. Using the first column as 'smiles'.")
        smiles_list = pd.read_csv(args.input_file, header=None).iloc[:, 0].tolist()
    except pd.errors.ParserError:
        print(f"Error: The input file is not a valid CSV")
        sys.exit(1)

    results = []
    for smiles in smiles_list:
        smiles = smiles.strip()
        if smiles:
            result = lipinski_pass(smiles)
            results.append({'smiles': smiles, 'Pass': result})

    results_df = pd.DataFrame(results)
    passed_smiles = results_df[results_df['Pass'] == True]['smiles']
    passed_smiles.to_csv(args.output, index=False, header=False)

if __name__ == '__main__':
    main()