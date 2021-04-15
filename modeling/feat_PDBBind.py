#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 15:13:38 2018

This script featurizes all PDBBIND complexes. Full length protein with 
partial charge calculated by AutodockTool is used (as will be in main training)
This script takes the cvset index (1-5) and the cvset (training, validation, test)
as input to create the .hdf files for the corresponding parition.

@author: cl321
"""
from tfbio.data import Featurizer
import pybel
import h5py
import pandas as pd
import numpy as np
import pickle
import os
import subprocess
import random
import shutil
import glob
import sys

config = sys.argv[1]
cvset = sys.argv[2]

biochem_df = pd.read_csv("../resources/PDBBIND_cvsplit1127.csv") 
print('total number of pairs', len(biochem_df))
datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/" 
wd = datasetdir + '/PDBBind2018/M2/I9/%s/'%config
if not os.path.exists(wd):
    os.mkdir(wd)
pdbpath = '/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/'

def unzip_lgd(lgdtmp):
    if not os.path.exists(lgdtmp+'/pose.df'):
        if not os.path.exists(lgdtmp):
            os.mkdir(lgdtmp)
        subprocess.call(['unzip',lgdtmp+'.zip','-d', lgdtmp])

def create_dataset(cvset): 
    # select the pairs in the parition
    cv_df = biochem_df.loc[biochem_df[config+'|cvset']==cvset]
    cv_df.reset_index(drop=True, inplace=True)
    outputpath = '{0}/{1}_set.hdf'.format(wd, cvset)

    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')
    
    with h5py.File(outputpath, 'w') as f:
        for idx, arow in cv_df.iterrows():
            
            geneid = arow['ENTREZ_GENE_ID']
            pdbid = arow['PDB_CODE']
            affinity = arow['RESULT_VALUE_min']
            
            if not os.path.exists('/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/PDBQT_noH/%s_PDBBIND.pdbqt'%pdbid): 
                print('path not found ', pdbid)
                continue
            
            # PDBBIND data is splitted between two folders so need to find the appropriate path
            if os.path.exists(pdbpath + '/v2018-other-PL/'+pdbid):
                pairpath = pdbpath + '/v2018-other-PL/'+pdbid +'/'
            elif os.path.exists(pdbpath + 'refined-set/'+pdbid):
                pairpath = pdbpath + '/refined-set/'+pdbid+'/'
            else:
                print('path not found ', pdbid)
                continue
       
            # read pocket into pybel 
            pocket = next(pybel.readfile('pdbqt', '/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/PDBQT_noH/%s_PDBBIND.pdbqt'%pdbid))
            pocket_coords, pocket_features = featurizer.get_features(pocket, 
                                                                     molcode=-1)
            assert (pocket_features[:, charge_idx] != 0).any() 
            
            # to get uniformedly processed partial charges, convert ligand mol2 to pdbqt using autodock 
            # for these pdbs ligand generated separated using fconv: 4kai, 3eyd, 4kb7, 4kbi
            # /n/groups/alquraishi/Apps/fconv_124_bin_linux32/fconv 4kai.pdb -l --t=4kai_fconv.mol2
            lgdpdbqt = pairpath + pdbid + '_ligand.pdbqt'
            lgdmol2 = pairpath + pdbid+'_ligand.mol2'
            # if no PDBQT, generate PDBQT
            if not os.path.exists(lgdpdbqt):               
                if 'mer' not in arow['ligand']:
                    cmd=["/home/cl321/bin/mglpython",
                        "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
                        "-l", lgdmol2, '-o', lgdpdbqt]
                else:
                    cmd=["/home/cl321/bin/mglpython",
                        "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py", 
                        "-r", lgdmol2, '-o', lgdpdbqt]                
                subprocess.call(cmd)      
            # read ligand into pybel     
            try:
                mol = next(pybel.readfile('pdbqt', lgdpdbqt))
            except IOError:
                with open('error_log.txt', 'a') as e:
                    e.write('file ' + pairpath + '\n')
                
            ligand_coords, ligand_features = featurizer.get_features(mol, 
                                                                     molcode=1)
            # originally checks for partial charge == 0, but it makes sense for benzene
            # to have to partial charge of 0 on the C (4w52) so instead I will
            # just keep a record of what the 0 partial charge molecules are
            try:
                assert (ligand_features[:, charge_idx] != 0).any()
            except AssertionError:
                with open('error_log.txt', 'a') as e:
                    e.write('partialcharge ' + pairpath + '\n')
            
            centroid = ligand_coords.mean(axis=0)
            ligand_coords -= centroid
            pocket_coords -= centroid
    
            data = np.concatenate((np.concatenate((ligand_coords, pocket_coords)),
                                np.concatenate((ligand_features, pocket_features))), axis=1)
            
            dataset = f.create_dataset(geneid+'|'+pdbid+'_xtal', data=data, shape=data.shape, 
                                       dtype='float32', compression='lzf')
            dataset.attrs['affinity'] = affinity

if __name__ == "__main__":
    create_dataset(cvset)
