#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:59:30 2019

This script selects the negative example poses for training (poses far from the
active site). I preselected a pose for each pair so that the dataset will be 
consisted of 50% poses far from the active site ('medium' in the table, <20A) 
and very far from the active site ('high' in the table, >20A). When such pose 
does not exist for a pair, the script selects what is available. 

@author: cl321
"""
import pandas as pd
import numpy as np
import sys
import random
import os
import subprocess
import shutil

begin = int(sys.argv[1])
end = int(sys.argv[2])
randomseed = 0
random.seed(randomseed)

datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/"    
wd = datasetdir + '/D13/M2/I2/'
expdir = datasetdir + '/symlink_Exp1/'  
# this is the table containing the pose annotations
biochem_df = pd.read_csv("../examples/DTC_cvsplit1127_dist.csv",  dtype={'ENTREZ_GENE_ID': object})
dist_dir = datasetdir + '/distance_from_ctr/'

pickdict = {'ENTREZ_GENE_ID':[], 'INCHIKEY_DESALT':[], 
            'max_modelname':[], 'max_poseID':[], 
            'distance':[], 'expt':[], 
            'dataset':[]}

ds_map = {'KD':'1', 'KI':'2'}

biochem_df_updated = biochem_df.copy()
for idx, arow in biochem_df.iloc[begin:end].iterrows():
    agene = arow['ENTREZ_GENE_ID']
    algd = arow['INCHIKEY_DESALT']
    ds = 'NEG' + ds_map[arow['dataset']]
    
    # read in the distance table for distance information
    pair_rmsd = pd.read_csv(dist_dir + agene + '/' + algd + '_distance.csv')
    
    # if the pair should have a 'medium' negative sample, sample a pose between
    # 10A and 20A and give an experimental value of log(2.7). If such pose doesn't 
    # exist for this pair, select a pose > 20A, but give it an expt value log(3.2).
    # record the final samples.
    if arow['dist_bin'] == 'medium':
        sel = pair_rmsd.loc[(10<pair_rmsd['distance'])&(pair_rmsd['distance']<20)]
        expt_aff = 2.7
        if len(sel) == 0:
            sel = pair_rmsd.loc[pair_rmsd['distance']>20]
            expt_aff = 3.2
            biochem_df_updated.loc[idx, 'dist_bin'] = 'high'
    # vice versa for the 'high' negative sample
    elif arow['dist_bin'] == 'high':
        sel = pair_rmsd.loc[pair_rmsd['distance']>20]
        expt_aff = 3.2
        if len(sel)==0:
            sel = pair_rmsd.loc[(10<pair_rmsd['distance'])&(pair_rmsd['distance']<20)]
            expt_aff = 2.7
            biochem_df_updated.loc[idx, 'dist_bin'] = 'medium'
    elif pd.isnull(arow['dist_bin']):
        continue
    
    if len(sel)>=1: 
        [rmsd, max_modelname, max_poseID] = sel.sample(1).iloc[0].tolist()
    # note error by giving NaNs to the pose
    else:
        print('no pose with chosen range rmsd ', agene, algd, arow['dist_bin'])
        rmsd = np.NaN
        max_modelname = np.NaN
        max_poseID = np.NaN
    
    pickdict['ENTREZ_GENE_ID'].append(agene)
    pickdict['INCHIKEY_DESALT'].append(algd)
    pickdict['max_modelname'].append(max_modelname)
    pickdict['max_poseID'].append(max_poseID)
    pickdict['distance'].append(rmsd)
    pickdict['expt'].append(expt_aff)
    pickdict['dataset'].append(ds)

# save both the pose and the modified table 
pickfile = pd.DataFrame(pickdict)
pickfile.to_csv(wd + '/picked_records_neg_only_{0}-{1}.csv'.format(begin, end), 
                index=False)
biochem_df_updated.iloc[begin:end].to_csv(wd+"DTC_cvsplit1127_dist_{0}-{1}.csv".format(begin, end),
                         index=False)
