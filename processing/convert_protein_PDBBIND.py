#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:19:25 2020

This script converts all proteins in PDBBIND2018 into pdbqt (H removed, 
partial charge calculated), which will be used for featurization
For parallelization, the conversion is performed in batches, so the script takes
the start and end index of a batch as input

@author: cl321
"""

import pandas as pd
import numpy as np
import pickle
import h5py
import glob
import pybel
from biopandas.pdb import PandasPdb
import subprocess
import os
from tfbio.data import Featurizer, make_grid, rotate
import sys

begin = int(sys.argv[1])
end = int(sys.argv[2])

# read and process PDBBIND pairs
general_raw = pd.read_table('INDEX_general_PL_data.2018', header=5)
general = general_raw['# =============================================================================='].str.split(expand=True)
general.rename(columns={0:'PDB_CODE', 1:'resolution', 2:'year', 
                          3:'-log(ki/kd)', 4:'affinity', 7:'ligand'}, inplace=True)
general.reset_index(drop=True, inplace=True)

# select the indicated batch
biochem = general.iloc[begin:end]

obabel = '/n/groups/alquraishi/Changchang/openbabel-2/bin/obabel'
pdbqt_dir = '/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/PDBQT_noH/'
                    
def desalt_PDB(inputfile, outputfile):
    with open(inputfile) as oldfile, open(outputfile, 'w') as newfile:
        for line in oldfile:
            if 'HETATM' not in line:
                newfile.write(line)

def autodock_pdbqt(inputfile, outputfile):
    pdbqt_cmd = ["/home/cl321/bin/mglpython",
                '/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py', 
                '-r', inputfile, 
                '-o', outputfile]
    subprocess.call(pdbqt_cmd)
    return outputfile

def extract_single_chain(xtal_path, pdbid):
    """sometimes the fulllength protein in PDBBIND contains multiple chains. 
    I will use the chain kept in the pocket file and extract that"""
    
    ppdb = PandasPdb()
    #identify the single to extract based on pocket.pdb
    pocket_df_pdbbind = ppdb.read_pdb(xtal_path.split('protein.pdb')[0]+'pocket.pdb').df
    pocket_chain = pocket_df_pdbbind['ATOM']['chain_id'].drop_duplicates().tolist()
    
    # select the single chain from full protein and save it
    fullprotein = ppdb.read_pdb(xtal_path).df
    protein_single_chain=PandasPdb()
    protein_single_chain.df['ATOM'] = fullprotein['ATOM'].loc[fullprotein['ATOM']['chain_id'].isin(pocket_chain)]
    chain_path = '{0}/{1}_chain{2}.pdb'.format(pdbqt_dir, pdbid, ','.join(pocket_chain))
    protein_single_chain.to_pdb(path=chain_path,records=None, gz=False, append_newline=False)
    return chain_path

def convert_PDBQT(PDBID):
    """main script to convert .pdb into .pdbqt"""
    try:
        xtal_path = glob.glob('/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/*/{0}/{0}_protein.pdb'.format(PDBID))[0]
    except IndexError:
        return 'noPDB'
    xtal_noH = pdbqt_dir + 'noH-PDBBIND_'+ PDBID+'.pdb'
    pdbqt_path = pdbqt_dir + PDBID+'_PDBBIND.pdbqt'
    desalted = pdbqt_dir + 'desalted-PDBBIND_'+PDBID+'.pdb'
    
    # first try remove H and convert to PDBQT
    removeH_cmd = [obabel, '-ipdb', xtal_path, '-opdb', '-O'+xtal_noH, '-d']
    subprocess.call(removeH_cmd)

    autodock_pdbqt(xtal_noH, pdbqt_path)
    
    if os.path.exists(pdbqt_path):
        return 'protein'
    else:
        # if the simple approach fails, desalt the protein first
        desalt_PDB(xtal_noH, desalted)
        autodock_pdbqt(desalted, pdbqt_path)   
        if os.path.exists(pdbqt_path):
            return 'desalted'
        else:
            # if desalting fails, extract single chain protein and repeat desalting and conversion
            single_chain_path = extract_single_chain(xtal_path, PDBID)
            removeH_cmd = [obabel, '-ipdb', single_chain_path, '-opdb', '-O'+xtal_noH, '-d']
            subprocess.call(removeH_cmd)
            desalt_PDB(xtal_noH, desalted)
            autodock_pdbqt(desalted, pdbqt_path)   
            if os.path.exists(pdbqt_path):
                return 'singlechain'
            else:
                return 'failed'

# perform conversion. Use the status returned as record 
biochem['status'] = biochem['PDB_CODE'].map(convert_PDBQT)
biochem.to_csv('PDBBIND_pdbqt_{0}-{1}.csv'.format(begin, end), index=False)
