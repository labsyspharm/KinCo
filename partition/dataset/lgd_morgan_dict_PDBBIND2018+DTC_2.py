#!/home/cl321/miniconda2/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 14:09:43 2018

This is the code to calculate pairwise distance between ligands. Because there 
are many ligands and the pairwise distance matrix will be huge, we will store
each row (the distance between a reference compound to all the compounds in the
dataset)in an individual file under that compound. 

@author: cl321
"""
import pandas as pd
#rdkit 2018.09.3
import rdkit
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import Chem

import sys
import pickle
import os
import numpy as np
import glob

biochem_df = pd.read_csv('../../resources/PDBBIND2018_DTC_1127.csv')
biochem_df.drop_duplicates(subset='INCHIKEY_DESALT', inplace=True)

wd = '.'
obabel  = '/n/groups/alquraishi/Changchang/openbabel-2a/bin/obabel'
import subprocess

sanitize_flags = Chem.rdmolops.SanitizeFlags.SANITIZE_ALL & ~Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES     

def ReadSmiles(smi):
    # first read in the molecule then apply sanitization flag to control what to sanitize for
    m = Chem.MolFromSmiles(smi, False)
    try:
        Chem.SanitizeMol(m, sanitize_flags)
    except ValueError:
        return None
    return m

def Import_FingerPrints(smi):
    m = ReadSmiles(smi)
    if m==None:
        m = smitocansmiles(smi)
    fp = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
    fp_vect = fp.ToBitString()
    return fp_vect

def smitocansmiles(smi):
    # if can't kekulize, convert explicitly to canonical smiles and read again
    tmp = wd + '/test/can.smi'
    cmd = [obabel, '-:'+smi, 
       '-ocan', '-O'+tmp]
    subprocess.call(cmd)
    smi_raw = pd.read_table(tmp, header = None)
    smi = smi_raw.iloc[0][0]
    
    m = Chem.MolFromSmiles(smi, False)
    Chem.SanitizeMol(m, sanitize_flags)
    return m

fp_dict = {}
print('generating fingerprint matrix...')
for idx, arow in biochem_df.iterrows():
    acmpd = arow['CANONICAL_SMILES_DESALT']
    inchikey = arow['INCHIKEY_DESALT']
    try:
        fp_dict[inchikey] = np.asarray(list(Import_FingerPrints(acmpd))).astype(int)
    except ValueError:
        print(arow['INCHIKEY_DESALT'], arow['CANONICAL_SMILES_DESALT'])
with open('morganfp_2_1024_dict_PDBBIND2018+DTC.p', 'wb') as handle:
    pickle.dump(fp_dict, handle)
print('saving fingerprint dictionary')
