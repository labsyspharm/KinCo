#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:33:53 2019

This script featurizes all the poses of each given kinase-compound pair into a 
hdf file. To parallelize, the dataset is batched. 

Input: 
    overwrite: if a .hdf file or energy predict exists, whether to regenerate features
    tgtgeneID_list: list of gene IDs to process
    lgdinchikey_list: list of compound inchikeys to process (must be in the same
                                                             order of tgtgeneID_list)
    machine: eg. exx2, under which GPU will this feature be predicted
    gpu: eg. 0, under which node will this feature be predicted
Ouput:
    hdf files with all poses of the kinase-compound pairs input

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
import sys
import gzip
import traceback
import logging
import time
import glob

overwrite = sys.argv[1]
tgtgeneID_list = sys.argv[2]
lgdinchikey_list = sys.argv[3]
machine = sys.argv[4]
gpu = sys.argv[5]

randomseed = 0
random.seed(randomseed)

datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/"    
expdir = datasetdir + '/symlink_Exp1/'
tmp_dir = '/n/scratch3/users/c/cl321/Datasets/D13/M2/I2/CNN_features/'
sfsdir = tmp_dir + machine + '/' + gpu +'/'
wd = datasetdir + '/D13/M2/I2/'
energy_path = wd + '/energy/'

geneidtokinasename=pickle.load(open(("../resources/GeneIDtoKinasename.p"),'rb'))

def unzip_lgd(lgdtmp):
    if os.path.exists(lgdtmp):
        shutil.rmtree(lgdtmp)
    if not os.path.exists(lgdtmp):
        os.mkdir(lgdtmp)
    subprocess.call(['unzip',lgdtmp+'.zip','-d', lgdtmp])

def read_mol(lgdtmp, lgdinchikey, tgtgeneID):
    # pdbqt format has to be used here in order to match the atom features with the coordinates of 
    # the docked structures 
    mol = next(pybel.readfile('pdbqt',  lgdtmp+'/PDBQT/'+lgdinchikey+'.pdbqt'))
    return mol

def read_pocket(tgtgeneID, modelname, subdomain):
    # read the homology model structure
    tgttmp = expdir+tgtgeneID     
    modelpath=(tgttmp +'/PDBQT/'+ modelname+'.pdbqt')
    pocket = next(pybel.readfile('pdbqt', modelpath))
    return pocket

def feat_pair(lgdinchikey, lgdchain, ligand_features, 
              tgtgeneID, subdomain, modelname, 
              featurizer, charge_idx):
    
    """
    This function combines the features of one homology model and one pose. It
    1) reads the homology model of the kinase; 
    2) loads the structure of the compound before docking; 
    3) replaces the ligand coordinates before docking with the coordinates 
    of the docked poses stored; 
    4) centers the complex around the pose; 
    5) concatenates the homology and the pose atoms
    """
    # 1) reads the homology model of the kinase
    pocket = read_pocket(tgtgeneID, modelname, subdomain)
    # 2) loads the structure of the compound before docking; 
    lgd_woH = lgdchain[~lgdchain['atom_type'].str.contains('H.')]
    
    try:
        pocket_coords, pocket_features = featurizer.get_features(pocket, molcode=-1)
        assert pocket_features.shape[0] == pocket_coords.shape[0], 'mismatched kinase features and coordinates'
        assert (pocket_features[:, charge_idx] != 0).any() 
        
        # 3) replaces the ligand coordinates before docking with the coordinates 
        ligand_coords = lgd_woH[['x', 'y','z']].as_matrix()
        assert ligand_features.shape[0] == ligand_coords.shape[0], 'mismatched ligand features and coordinates'
        assert (ligand_features[:, charge_idx] != 0).any()
        
    except AssertionError:
        print ('mismatched features and coordinates ', tgtgeneID, lgdinchikey)
        return None
    
    # 4) centers the complex around the pose; 
    centroid = ligand_coords.mean(axis=0)
    ligand_coords -= centroid
    pocket_coords -= centroid
    
    # 5) concatenates the homology and the pose atoms
    coords = np.concatenate((ligand_coords, pocket_coords))
    features = np.concatenate((ligand_features, pocket_features))
    
    return [coords, features]
    
def df_to_hdf(tgtgeneID, lgdinchikey):
    """
    This is the main function to featurize a kinase-compound pair
    
    """
    # set up directory 
    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')
    
    tgttmp = expdir+tgtgeneID
    lgdtmp = tgttmp + '/' +lgdinchikey
    pairbase = lgdtmp + '/'
    unzip_lgd(lgdtmp)
    
    sfstgtdir = sfsdir + tgtgeneID + '/'
    if not os.path.exists(sfstgtdir):
        os.mkdir(sfstgtdir)
    
    # first generate prime.hdf files so that we can differentiate the unsucessful
    # .hdfs from those that were finished and only transfer the latter
    outputpath = sfstgtdir + lgdinchikey+'_prime.hdf'
    
    all_coords = []
    all_features = []
    all_names = []
    
    with h5py.File(outputpath, 'w') as f:
        try:
            subdomain = geneidtokinasename[tgtgeneID] + '_D0'
        except IOError:
            print ('wrong key ', tgtgeneID)
        
        # load poses
        pose_df = pd.read_pickle(pairbase + 'pose.df')
        pose_df = pose_df.drop_duplicates(subset=['ModelName', 'PoseID'])
        pose_df = pose_df.dropna(subset=['ModelName', 'PoseID'], axis = 0)

        # reference ligand file is the same for all poses, so just load it 
        # once upfront
        try:
            mol = read_mol(lgdtmp, lgdinchikey, tgtgeneID)
            _, ligand_features = featurizer.get_features(mol, molcode=1)
        except IOError:
            print('no ligand pdbqt')
            return None
        
        # vary the coordinates for the different poses and featurize them
        for index, arow in pose_df.iterrows():
            ligandchain = arow['LigandChain']
            modelname = arow['ModelName']
            
            try:
                [pair_coords, pair_features] = feat_pair(lgdinchikey, ligandchain, ligand_features,
                          tgtgeneID, subdomain, modelname, 
                          featurizer, charge_idx)
            except IOError:
                print('no target pdb')
                return None
            
            data = np.concatenate((pair_coords, pair_features), axis=1)
            name = modelname + '|' + arow['PoseID']

            f.create_dataset(name, data=data, shape=data.shape, dtype='float32', compression='lzf')     
    
    # makesure file is generated and clean up
    if os.path.exists(outputpath):
        size = os.path.getsize(outputpath)/1000./1000
        if size == 0:
            print('ERROR: hdf file not produced', tgtgeneID, lgdinchikey)
        else:
            shutil.move(outputpath, sfstgtdir + lgdinchikey+'.hdf')
            
    if (os.path.exists(lgdtmp+'.zip')) and (os.path.getsize(lgdtmp+'.zip') > 2000):
        shutil.rmtree(lgdtmp)
        
    return [all_coords, all_features, all_names]

if __name__ == "__main__":
    
    tgtgeneIDlist = tgtgeneID_list.split(',')
    lgdinchikeylist = lgdinchikey_list.split(',')
    
    assert len(tgtgeneIDlist) == len(lgdinchikeylist)
    
    # for each kinase-compound pair, check if it has been predicted or featurize
    # featurize if the pair has not been processed
    for atgtgeneID, algdinchikey in zip(tgtgeneIDlist, lgdinchikeylist):
        print ('working on ', atgtgeneID, algdinchikey)
        
        energylgd = glob.glob(energy_path + '/*/'+atgtgeneID + '/' + algdinchikey+'_tf.csv')
        sfslgdhdf = glob.glob('{0}/CNN_features/*/*/{1}/{2}.hdf'.format(tmp_dir, atgtgeneID, algdinchikey))
        check_files = (len(sfslgdhdf)==0) and (len(energylgd)<5)
        if  check_files or (overwrite == 'True'):
            df_to_hdf(atgtgeneID, algdinchikey)
            
