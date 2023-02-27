#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 18:07:32 2019

This script calculates RMSD to the center for all poses of a kinase-compound pair
The dataset is batched so the script takes as input the beginning and end index
of the batch and iterates through each of the pairs

@author: cl321
"""
import subprocess 
import os
import numpy as np
import pandas as pd
import shutil
import sys

DeepAlign='/home/cl321/DeepAlign/DeepAlign'
parentfldr = '/n/groups/alquraishi/Changchang/Datasets/symlink_Exp1/'
distance_dir = '/n/groups/alquraishi/Changchang/Datasets/distance_from_ctr/'

begin = int(sys.argv[1])
end = int(sys.argv[2])

def unzip_lgd(lgdtmp):
    if not os.path.exists(lgdtmp+'/pose.df'):
        if not os.path.exists(lgdtmp):
            os.mkdir(lgdtmp)
        status = subprocess.call(['unzip',lgdtmp+'.zip','-d', lgdtmp])
        if status == 9:
            raise FileNotFoundError
    return lgdtmp+'/'

def get_ref_coordinates(lgd_inchikey, tgt_geneid, hm):
    """ 
    This function determines the center of the active site of a kinase structure
    (homology model) by aligning its structure to a reference kinase structure 
    and using the residues in the active site of the reference structure to 
    infer the coordinates of the center of the query kinase structure
    
    input:
        lgd_inchikey: inchikey of the compound
        tgt_geneid: gene ID of the kinase of interest
        hm: name of the homology model structure
    output:
        coordinates of the center of the active site of the query kinase structure
    """
    
    tgtfldr = parentfldr + '/' + tgt_geneid + '/'
    lgdfldr = tgtfldr + '/' + lgd_inchikey + '/'
    
    os.chdir(lgdfldr)
    
    # align query homology model structure to the reference structure
    # find the reference target file 1a9u .1
    siteref=('../resources/1a9u.1.1.1.tgt.pdb')
    recordaln= tgt_geneid + '-' + lgd_inchikey + '-' + hm +'_align'
    outputaln = recordaln+'.pdb'
    parsealign= lgdfldr + hm + '-align_parsed.pdb'
    tgt_location = tgtfldr +'/PDBQT/'+ hm +'.pdbqt'
    
    if not os.path.isfile(parsealign):
        #runs DeepAlign to superimpose 1A9U onto query structure
        subprocess.call([DeepAlign, (siteref), (tgt_location), "-o", (recordaln)])
        #.pdb file actually contains both the pdb of the reference and model, so 
        #here I extract only the structure of the model and put it under its target
        #file
        with open((outputaln)) as f:
            lines=f.readlines()
            with open(parsealign, 'a+') as g:
                for line in lines:
                    if 'TER' in line:
                        break
                    g.write(line)
        #delete the temporary files
        outputalnlist=[(recordaln+'.fasta'), (recordaln+'.local'),
                       (recordaln+'.score'), (recordaln+'.pdb')]
        for anoutput in outputalnlist:    
            os.remove(anoutput)
    
    # define binding site using R515 and M865 in 1A9U
    # find the coordinates of these two residues in the superimposed 1A9U structure
    # and calculate the mean of the new coordinates as the center coordinates
    with open(parsealign,'r') as g:
        lines=g.readlines()
        for line in lines:
            alist=line.split()
            if alist[1]=="865" and alist[2]=="O" and alist[3]=="MET" and alist[4]=="A":
                try:
                    coord1 = np.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord1=np.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord1=np.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break
        for line in lines:
            alist=line.split()
            if alist[1]=="515" and alist[2]=="NH2" and alist[3]=="ARG" and alist[4]=="A":
                try:
                    coord2 = np.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord2=np.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord2=np.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break        
    center = (coord1+coord2)/2
    
    return center

def calc_RMSD(refpt, querydf):
    """
    calculates the distance of the query ligand to the center point the
    distance between two points (represents the ligand as its centroid and
    calculate the distance between the centroid and the center of active site)
    
    input: 
        refpt: [x, y, z] coordinates of the active site center
        querydf: dataframe with the coordinates of the ligand
    output:
        the distance between the two points
    """
    query_p = querydf.loc[~querydf['atom_type'].str.contains('H')]
    querycoord = query_p[['x','y','z']].values
    query_pt = querycoord.mean(axis=0)
    rmsd = np.sqrt(sum((query_pt - refpt)**2))
    return rmsd

biochem_df = pd.read_csv('../resources/DTC_cvsplit1127.csv', dtype={'ENTREZ_GENE_ID': object})
biochem_df.drop_duplicates(subset = ['ENTREZ_GENE_ID', 'INCHIKEY_DESALT'], inplace=True)
biochem_df.reset_index(drop=True, inplace=True)

for idx, apair in biochem_df.iloc[begin:end].iterrows(): 
    tgt_geneid = apair['ENTREZ_GENE_ID']
    lgd_inchikey = apair['INCHIKEY_DESALT']
    
    # set up the file directory. skip the pair if distance file already exists
    if not os.path.exists(distance_dir+tgt_geneid+'/'):
        os.mkdir(distance_dir+tgt_geneid+'/')
    
    if os.path.exists(distance_dir+tgt_geneid+'/'+lgd_inchikey+'_distance.csv'):
        continue
    
    pairpath = parentfldr + tgt_geneid+'/'+lgd_inchikey
    
    # load the coordinates of all poses
    try:
        pairdir = unzip_lgd(pairpath)
        posedf = pd.read_pickle(pairdir+'pose.df')
    except FileNotFoundError:
        with open(distance_dir+ 'error_file.txt', 'a+') as ef:
            ef.write(tgt_geneid+' '+lgd_inchikey+'\n')
            continue
        
    # There are multiple poses for one homology model. To save compute, calculate
    # the center of the homology models first and refer to them later on (instead
    # of repeating the same calculation for each pose docked into the same homology
    # model)
    hm_center_dict = {}
    for ahm in posedf['ModelName'].drop_duplicates().tolist():
        aln_center = get_ref_coordinates(lgd_inchikey, tgt_geneid, ahm)
        hm_center_dict[ahm] = aln_center
    
    # calculate the distance of each pose to the center of the active site
    # as distance between two points
    res_dict = {'modelname':[], 'poseID':[], 'distance':[]}
    
    for idx, arow in posedf.iterrows():
        modelname = arow['ModelName']
        poseID = arow['PoseID']
        ligand_chain = arow['LigandChain']
        pose_RMSD = calc_RMSD(hm_center_dict[modelname], ligand_chain)
        
        res_dict['modelname'].append(modelname)
        res_dict['poseID'].append(poseID)
        res_dict['distance'].append(round(pose_RMSD, 3))
        
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(distance_dir+tgt_geneid+'/'+lgd_inchikey+'_distance.csv', 
                  index=False)
    
    shutil.rmtree(pairdir)
