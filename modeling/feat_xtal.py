#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 15:13:38 2018

This script featurizes all the kinase-compound pairs in PDBBIND2018. Note that
the kinase structure were processed to remove H and partial charge calculated
using Autodock tool scirpt by converting to PDBQT. The full length protein is 
used, instead of just the pocket as in Pafnucy. 

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

biochem_df_all = pd.read_csv("/n/groups/alquraishi/Changchang/Datasets/D13/dataset/PDBBIND-kinase_DTC_cvsplit1127.csv",
                         dtype={'ENTREZ_GENE_ID':object}) 
biochem_df = biochem_df_all.loc[biochem_df_all['dataset']=='PDBBIND']

datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/" 
wd = datasetdir + '/D13/M2/I1/'
# use processed PDBQT files
pdbqt_path = '/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/PDBQT_noH/' 
pdbpath = '/n/groups/alquraishi/Changchang/Datasets/PDBBind2018/'
configs = ['cv1', 'cv2', 'cv3', 'cv4', 'cv5']

for i in configs:
    if not os.path.exists(i):
        os.mkdir(i)

def create_dataset(): 
        
    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')
    
    
    for idx, arow in biochem_df.iterrows():
            
        geneid = arow['ENTREZ_GENE_ID']
        pdbid = arow['PDB_CODE']
        affinity = arow['RESULT_VALUE_min']

        # use PDBQT with H removed
        pocketpdbqt = pdbqt_path + pdbid+'_PDBBIND.pdbqt'       
        # read pocket into pybel 
        pocket = next(pybel.readfile('pdbqt', pocketpdbqt))
        
        pocket_coords, pocket_features = featurizer.get_features(pocket, 
                                                                 molcode=-1)
        assert (pocket_features[:, charge_idx] != 0).any() 
        
        if os.path.exists(pdbpath + '/v2018-other-PL/'+pdbid):
            pairpath = pdbpath + '/v2018-other-PL/'+pdbid +'/'
        elif os.path.exists(pdbpath + 'refined-set/'+pdbid):
            pairpath = pdbpath + '/refined-set/'+pdbid+'/'
        else:
            print('path not found ', pdbid)
            continue
        
        # to get uniformedly processed partial charges, convert ligand mol2 to pdbqt using autodock 
        # for these pdbs ligand generated separated using fconv: 4kai, 3eyd, 4kb7, 4kbi
        # /n/groups/alquraishi/Apps/fconv_124_bin_linux32/fconv 4kai.pdb -l --t=4kai_fconv.mol2
        lgdpdbqt = pairpath + pdbid + '_ligand.pdbqt'
        lgdmol2 = pairpath + pdbid+'_ligand.mol2'
        # if no PDBQT, generate PDBQT
        if not os.path.exists(lgdpdbqt):
            # most ligands are small molecules
            if 'mer' not in arow['ligand']:
                cmd=["/home/cl321/bin/mglpython",
                    "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
                    "-l", lgdmol2, '-o', lgdpdbqt]
            # some ligands are peptides
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
                        
        for acv in configs:
            cvset = arow['%s|cvset'%acv]
            
            if cvset in ['training', 'validation', 'test']:
                cvpath = '{0}/{1}/{2}_set_xtal.hdf'.format(wd, acv, cvset)
            elif pd.isnull(cvset):
                continue
            else:
                print('cvset not recognized')
                continue
            
            with h5py.File(cvpath, 'a') as f:
                dataset = f.create_dataset(geneid+'|'+pdbid+'_PDBBIND', data=data, shape=data.shape, 
                                           dtype='float32', compression='lzf')
                dataset.attrs['affinity'] = affinity

if __name__ == "__main__":
    create_dataset()
