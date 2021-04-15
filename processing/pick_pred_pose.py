#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:59:30 2019

This script picks out the pose with the lowest energy considered by both vina 
and the energy model (lowest average rank of vina and M1) and records it in a 
file. We assume that the compounds are ATP-competitive inhibitors so we will
focus on only the poses within 10A of the active site. 

For parallelization, the dataset is batched, so this script takes the beginning
and the end indexes of the batch. It also requires specifying which models to 
use (energyprefix). The script uses the energy prediction file from the model (M1),
the prediction file by vina, and the file with distance from active site for each pose.

@author: cl321
"""
import pandas as pd
import numpy as np
import sys
import random

randomseed = 0
random.seed(randomseed)

datasetdir =  "/n/groups/alquraishi/Changchang/Datasets/"    
wd = datasetdir + '/D13/M2/I2/'
expdir = datasetdir + '/symlink_Exp1/'  
biochem_df = pd.read_csv("../resources/DTC_cvsplit1127.csv",  dtype={'ENTREZ_GENE_ID': object})
dist_dir = datasetdir + '/distance_from_ctr/'

begin = sys.argv[1]
end = sys.argv[2]
energyprefix = sys.argv[3]

sfsenergy = '{0}/energy/{1}/'.format(wd, energyprefix)
vinadir = datasetdir + 'vina_energy/'

def Get_MinAff(pred_df_vina, pred_df, rmsd_df):
    """
    To select the pose with the lowest average rank between vina and the model
    by ranking the predictions from each model and taking the average 
    
    input:
        pred_df_vina: prediction table by vina containing all poses 
        pred_df: prediction table by model (e.g. M1) containing all poses 
        rmsd_df: table containing the RMSD of all poses from the active site
    output: 
        a list. Information on the pose with the lowest average rank incluidng
        homology model name, poseID, the average rank, prediction from model,
        prediction from vina and the distance from the active site
    """
    pred_df.drop_duplicates(inplace=True)
    rmsd_df.drop_duplicates(inplace=True)
    pred_df_vina.drop_duplicates(inplace=True)

    # merge the three tables (model, vina, and RMSD) into one 
    pred_rmsd = pred_df.merge(rmsd_df)
    pred_vina_rmsd = pred_rmsd.merge(pred_df_vina)

    # get the rank of each pose based on vina and model. take the average
    pred_vina_rmsd['rank_vina'] = pred_vina_rmsd['energy_kcal'].rank(method = 'min')
    pred_vina_rmsd['rank_model'] =  pred_vina_rmsd['Kdmodel'].rank(method = 'min')
    pred_vina_rmsd['rank_avg']  = pred_vina_rmsd[['rank_vina', 'rank_model']].mean(axis=1)
    # select the one with lowest average rank
    pred_vina_rmsd_srt = pred_vina_rmsd.sort_values(by='rank_avg').reset_index(drop=True)
    return pred_vina_rmsd_srt.iloc[0][['modelname', 'poseID', 'rank_avg', 'Kdmodel', 'energy_kd','distance']].tolist()

def pick_by_rank(pred_df_vina, pred_df_model, pred_rmsd):
    """
    To select the pose within the active site AND with the lowest average rank 
    between vina and the model by ranking the predictions from each model and 
    taking the average 
    
    input:
        pred_df_vina: prediction table by vina containing all poses 
        pred_df: prediction table by model (e.g. M1) containing all poses 
        rmsd_df: table containing the RMSD of all poses from the active site
    output: 
        a list. Information on the pose with the lowest average rank incluidng
        homology model name, poseID, the average rank, prediction from model,
        prediction from vina and the distance from the active site
    """
    pred_df_vina.drop_duplicates(inplace=True)
    pred_rmsd.drop_duplicates(inplace=True)
    pred_df_model.drop_duplicates(inplace=True)
    
    # merge the three tables (model, vina, and RMSD) into one 
    pred_model_rmsd = pred_df_model.merge(pred_rmsd)
    pred_model_rmsd_vina = pred_model_rmsd.merge(pred_df_vina)
    # select the poses within 10A of active site
    pred_activesite = pred_model_rmsd_vina.loc[pred_model_rmsd_vina['distance']<=10]
    
    # get the rank of each pose based on vina and model. take the average
    pred_activesite['rank_vina'] = pred_activesite['energy_kcal'].rank(method = 'min')
    pred_activesite['rank_model'] = pred_activesite['Kdmodel'].rank(method = 'min')
   
    # select the one with lowest average rank
    pred_activesite['rank_avg'] = pred_activesite[['rank_vina', 'rank_model']].mean(axis=1)
    pred_merge_rmsd_srt = pred_activesite.sort_values(by='rank_avg').reset_index(drop=True)
    return pred_merge_rmsd_srt.iloc[0][['modelname', 'poseID', 'rank_avg', 'Kdmodel', 'energy_kd', 'distance']]
    
pickdict = {'ENTREZ_GENE_ID':[], 'INCHIKEY_DESALT':[], 
            'max_modelname':[], 'max_poseID':[], 'kdmodel':[], 'kdvina':[],
            'prediction':[], 'expt':[],'distance':[], 
            'dataset':[]}

for idx, arow in biochem_df.iloc[int(begin):int(end)].iterrows():
    
    # initialize each row
    agene = arow['ENTREZ_GENE_ID']
    algd = arow['INCHIKEY_DESALT']
    ds = arow['dataset']
    expt = arow['RESULT_VALUE_min']
    
    energyfilepath = sfsenergy+agene+'/'+algd+'_tf.csv'
    vinafilepath = vinadir+agene+'/'+algd+'_energy.csv'
    
    # read model, vina, and RMSD tables for each pair 
    try:
        resultmodel = pd.read_csv(energyfilepath)
        resultvina = pd.read_csv(vinafilepath)
        pair_rmsd = pd.read_csv(dist_dir + agene + '/' + algd + '_distance.csv')
    except pd.errors.EmptyDataError:
        print(agene, algd)
        continue
    except IOError:
        print(agene, algd)
        continue        
    
    try:
        # find the pose with the lowest average rank and within 10A
        [min_modelname, min_poseID, min_pred_rank, minkdmodel, minkdvina, distance] = pick_by_rank(resultvina, resultmodel, pair_rmsd)
    # for pairs where no poses within 10A exists (only 1 pair), select the pose with the lowest average rank
    except ValueError:
        [min_modelname, min_poseID, min_pred_rank, minkdmodel, minkdvina, distance] = Get_MinAff(resultvina, resultmodel, pair_rmsd)
    except IndexError:
        [min_modelname, min_poseID, min_pred_rank, minkdmodel, minkdvina, distance] = Get_MinAff(resultvina, resultmodel, pair_rmsd)

    pickdict['ENTREZ_GENE_ID'].append(agene)
    pickdict['INCHIKEY_DESALT'].append(algd)
    pickdict['max_modelname'].append(min_modelname)
    pickdict['max_poseID'].append(min_poseID)
    pickdict['distance'].append(distance)
    pickdict['prediction'].append(min_pred_rank)
    pickdict['kdmodel'].append(minkdmodel)
    pickdict['kdvina'].append(minkdvina)
    pickdict['expt'].append(expt)
    pickdict['dataset'].append(ds)

pickfile = pd.DataFrame(pickdict)
pickfile.to_csv('{0}/picked_records_{1}-{2}.csv'.format(sfsenergy, begin, end), 
                index=False)
