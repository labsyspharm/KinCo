#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 15:13:38 2018
Featurizes the selected pose for each kinase-compound pair. For parallelization,
the kinase-compound pairs are batched and this script featurizes one batch. 
Therefore it requires the start and end indexes of the list to be processed.
For efficiency, featurizes for all 5 CVsets at the same time
Featurization script is based on Pafnucy prepare.py

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

begin = sys.argv[1]
end = sys.argv[2]

geneidtokinasename=pickle.load(open(("../resources/GeneIDtoKinasename.p"),'rb'))

datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/"    
ensemblerdir = "/n/groups/alquraishi/Changchang/Ensembler/" 
wd = datasetdir + '/D13/M2/I2/' # path of the working directory for each iteration
expdir = datasetdir + '/symlink_Exp1/'

configs = ['cv1', 'cv2', 'cv3', 'cv4', 'cv5']
for aconf in configs:
    if not os.path.exists(wd+aconf):
        os.mkdir(wd+aconf)
        
# to pick poses using the 5 models separately for training
biochem_df = pd.read_csv(wd+'picked_records_wneg.csv',
                         dtype={'ENTREZ_GENE_ID':object})

def unzip_lgd(lgdtmp):
    if not os.path.exists(lgdtmp+'/pose.df'):
        if not os.path.exists(lgdtmp):
            os.mkdir(lgdtmp)
        subprocess.call(['unzip',lgdtmp+'.zip','-d', lgdtmp])

def extract_feat(tgtgeneID,lgdinchikey, max_modelname, pose_df, max_poseID):
    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')
    
    # 1) read ligand into pybel 
    lgdtmp = expdir+tgtgeneID+'/'+lgdinchikey
    lgdpdbqt = lgdtmp+'/PDBQT/'+lgdinchikey+'.pdbqt'
    mol = next(pybel.readfile('pdbqt',  lgdpdbqt))
    
    # 2) read in target pocket
    # to get uniformedly processed partial charges, convert pocket to pdbqt using autodock
    tgtpdb = expdir+tgtgeneID+'/PDB/'+max_modelname+'.pdb'
    tgtpdbqt = expdir+tgtgeneID+'/PDBQT/'+max_modelname+'.pdbqt'
    
    # read target into pybel 
    pocket = next(pybel.readfile('pdbqt', tgtpdbqt))
    
    max_ligandchain = pose_df.loc[(pose_df['ModelName']==max_modelname) & 
                            (pose_df['PoseID']==max_poseID)].iloc[0]['LigandChain']
    max_lgd_woH = max_ligandchain[~max_ligandchain['atom_type'].str.contains('H.')]

    _, ligand_features = featurizer.get_features(mol, molcode=1)
    ligand_coords = max_lgd_woH[['x', 'y','z']].as_matrix()
    
    assert ligand_features.shape[0] == ligand_coords.shape[0], 'mismatched features and coordinates'
    assert (ligand_features[:, charge_idx] != 0).any()
    pocket_coords, pocket_features = featurizer.get_features(pocket, molcode=-1)
    assert (pocket_features[:, charge_idx] != 0).any() 
    
    # 3) center the complex around the ligand and combine ligand and target
    centroid = ligand_coords.mean(axis=0)
    ligand_coords -= centroid
    pocket_coords -= centroid
    
    data = np.concatenate((np.concatenate((ligand_coords, pocket_coords)),
                        np.concatenate((ligand_features, pocket_features))), axis=1)
    return data
                
def create_dataset(begin, end): 
    sel_biochem = biochem_df.iloc[int(begin):int(end)]
    # if a pair has both Kd, Kd measurements, want to retain both measurements 
    # as separate data points
    picked_group = sel_biochem.groupby(['ENTREZ_GENE_ID', 'INCHIKEY_DESALT'])
    
    for g_idx, agroup in picked_group:
        
        tgtgeneID = agroup.iloc[0]['ENTREZ_GENE_ID']
        lgdinchikey = agroup.iloc[0]['INCHIKEY_DESALT']
        
        if not os.path.exists(expdir+tgtgeneID):
            print('no tgt')
            continue
        
        lgdtmp = expdir+tgtgeneID+'/'+lgdinchikey
        pairbase = lgdtmp + '/'
    
        if not os.path.exists(lgdtmp):
            unzip_lgd(lgdtmp)
        
        lgdpdbqt = lgdtmp+'/PDBQT/'+lgdinchikey+'.pdbqt'
        
        pose_df = pd.read_pickle(pairbase + 'pose.df')
        
        # use both pdbqt for ligand and target to get partial charge
        if not os.path.exists(lgdpdbqt):
            
            print('generating lgd pdbqt', lgdtmp)
            
            # if no mol2 generate mol2 from SMILES
            if not os.path.isfile((lgdtmp +'/mol2/'+lgdinchikey +'.mol2')):
                if not os.path.exists(lgdtmp +'/mol2/'):
                    os.mkdir(lgdtmp +'/mol2/')
                with open((lgdtmp+'/smiles'), 'r+') as smi:
                    ligand_smiles=smi.readline().strip('\n')
                obabel =  '/n/groups/alquraishi/Changchang/openbabel-2a/bin/obabel'
                subprocess.call([obabel, ('-:'+ligand_smiles), '-omol2', '-O', 
                                    (lgdtmp +'/mol2/' + lgdinchikey + '.mol2'), '--gen3d'])
            
            # generate PDBQT
            if not os.path.exists(lgdtmp+'/PDBQT/'):
                os.mkdir(lgdtmp+'/PDBQT/')
            cmd=["/home/cl321/bin/mglpython",
            "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
            "-l", (lgdtmp +'/mol2/' + lgdinchikey +'.mol2'), '-o', lgdpdbqt]
            subprocess.call(cmd)
        
        # main featurization
        for index, arow in agroup.iterrows():
            
            try:
                affinity = arow['expt']
                data_type = arow['dataset']
                name = tgtgeneID + '|' + lgdinchikey + '_' + data_type
            
                for acv in configs:
                    max_modelname = arow['%s|max_modelname'%acv]
                    max_poseID = arow['%s|max_poseID'%acv]
                    cvset = arow['%s|cvset'%acv]
                    if cvset in ['training', 'validation', 'test']:
                        cvpath = '{0}/{1}/{2}_set_{3}-{4}.hdf'.format(wd, acv, cvset, begin, end)
                    elif pd.isnull(cvset):
                        continue
                    else:
                        print('cvset not recognized')
                        continue
                    
                    try:
                        with h5py.File(cvpath, 'a') as f:
                            data = extract_feat(tgtgeneID, lgdinchikey, max_modelname, pose_df, max_poseID)
                            dataset = f.create_dataset(name, data=data, shape=data.shape, dtype='float32', compression='lzf')
                            dataset.attrs['affinity'] = affinity
                    except AssertionError:
                        print('mismatched features and coordinates ', lgdtmp, acv)
                    except Exception:
                        logging.error(traceback.format_exc())
                        with open('{0}/{1}/exceptions.txt'.format(wd, acv), 'a') as e:
                            np.savetxt(e, [str(index), tgtgeneID, lgdinchikey], fmt='%s')
                        continue                        
        
            except Exception:
                logging.error(traceback.format_exc())
                with open(wd + '/exceptions.txt', 'a') as e:
                    np.savetxt(e, [str(index), tgtgeneID, lgdinchikey], fmt='%s')
                continue

        if os.path.exists(lgdtmp+'.zip'):
            shutil.rmtree(lgdtmp)    
                                    
if __name__ == "__main__":
    create_dataset(begin, end)
