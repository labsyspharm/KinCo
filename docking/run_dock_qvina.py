#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:58:18 2017

@author: cl321
"""

""" 
This script identifies the ligand-target pairs that need to be docked and 
submit homology model docking jobs for each of the pairs. The pairs to dock are
in the input file. 
 """

import pickle
import argparse
import os
import glob
import shutil
import numpy as np
import subprocess
import pandas as pd
import csv
import random 

np.random.seed(1234)

parser = argparse.ArgumentParser()
parser.add_argument('start', help='index to start')
parser.add_argument('end', help='index to stop')
parser.add_argument('ncpus', help='number of cpus')
args = parser.parse_args()

"""set up the directory structure"""
# the file where docking is performed
parentfldr = '/n/groups/alquraishi/Changchang/Ensembler/symlink_Exp1/' 
# the tmp dir to store intermdiate docking results
parenttmpfldr = '/n/scratch3/users/c/cl321/Ensembler/symlink_Exp1/'
# working directory with the scripts
wd = '/n/groups/alquraishi/Changchang/Ensembler/Exp31/'
scratch = '/n/scratch3/users/c//cl321/Ensembler/Exp31/'
# ultimate storage of the docked results
dataset =  '/n/groups/alquraishi/Changchang/Datasets/symlink_Exp1/'

"""loading"""
data = pd.read_csv('../resources/DTC_1127.csv', 
                       low_memory = False)
data.drop_duplicates(subset=['ENTREZ_GENE_ID','INCHIKEY_DESALT'], inplace=True)
print(len(data))
computationcost = 560
geneidtokinasename=pickle.load(open(("../resources/GeneIDtoKinasename.p"),'rb'))

def Setup_Fldrs(geneID, ligand_inchikey, ligand_smiles):
    """ 
    set up subdirectionary structure for a unique kinase compound pair 
    input: the gene ID, ligand inchikey and ligand SMILES string 
    (CANONICAL SMILES from the input file) 
    output: the path of the kinase-compound directory under tmp directory
   """
    # create folders for targets under the parent 
    # create repertoire target folder in docking dir so that we can store the homology
    # models for the target
    tgtfldr = parentfldr + '/' + geneID + '/'
    if not os.path.exists(tgtfldr):
        os.mkdir(tgtfldr)
    if not os.path.exists((tgtfldr+'/PDB')): 
        os.mkdir((tgtfldr+'/PDB'))
    if not os.path.exists((tgtfldr+'/PDBQT')):
        os.mkdir((tgtfldr+'/PDBQT'))
    
    tgttmpfldr = parenttmpfldr + '/' + geneID + '/'
    if not os.path.exists(tgttmpfldr):
        os.mkdir(tgttmpfldr)    
    
    # create ligand files under targets
    # initialize temp ligand fldr under targets
    lgdfldr = tgttmpfldr + '/'+ ligand_inchikey + '/'
    if not os.path.exists(lgdfldr+'smiles'):
        if not os.path.exists(lgdfldr):
            os.mkdir(lgdfldr)
        with open((lgdfldr+'smiles'), 'w+') as smi:
            smi.write(ligand_smiles)
    else:
        with open((lgdfldr+'smiles'), 'r+') as smi:
            line=smi.readline().strip('\n')
            smi.close()
        if line != ligand_smiles:
            with open((lgdfldr+'smiles'), 'w+') as smi:
                smi.write(ligand_smiles)
                smi.close()
    if not os.path.exists((lgdfldr+'/mol2')):
        os.mkdir((lgdfldr+'/mol2'))
    if not os.path.exists((lgdfldr+'/PDBQT')):
        os.mkdir((lgdfldr+'/PDBQT'))     
    if not os.path.exists(lgdfldr+'/poses/'):
        os.mkdir(lgdfldr+'/poses/')              
    return lgdfldr

def Load_RunNumber(lgdtmp, subdomain):
    """ 
    determine the number of dockings to be performed for each homology model
    To ensure kinases with different number of homology models available undergo 
    the same number of docking, the number of runs is calculated by dividing 
    the total computational cost by the number of homology models. When there 
    are more homology models available than computational cost, models are 
    randomly sampled
    input: the path of the kinase-compound directory and name of the kinase domain
    output: minimal number of runs for each homology model
    """
    modellist = runs_dict[subdomain]
    
    # calculate runs if have not
    if glob.glob(lgdtmp+'run_number=*') == []:
        if len(modellist) < computationcost:
                runs= np.floor(computationcost/float(len(modellist)))
                left = computationcost % len(modellist)
                overmodellist = np.random.choice(modellist, left, replace = False)
                undermodellist = [x for x in modellist if x not in overmodellist]
                np.savetxt((lgdtmp+'run_number='+str(int(runs))), undermodellist, fmt = "%s")
                np.savetxt((lgdtmp+'run_number='+str(int(runs+1))), overmodellist, fmt = "%s")
        else:
            runs = 1
            undermodellist = np.random.choice(modellist, computationcost, replace=False)
            np.savetxt((lgdtmp+'run_number='+str(int(runs))), undermodellist, fmt = "%s")
            overmodellist = []
            np.savetxt((lgdtmp+'run_number='+str(int(runs+1))), overmodellist, fmt = "%s")
    # load the number of runs if already calculated
    else: 
        tmp_runlist=[int(i.split('=')[-1]) for i in sorted(glob.glob(lgdtmp+'run_number=*'))]
        runs = min(tmp_runlist)
    return runs   
     
def Submit_Docking(geneID, ligand_inchikey, subdomain, runs):
    """ 
    submit a docking job for each kinase-compound pair. Command is printed
    and to be output to a command submission .txt on terminal so multiple pairs 
    can be submitted then
    input: kinase gene ID, ligand inchikey, Uniprot name of kinase domain, 
          number of times to perform docking
    output: None, but command will be printed
    """
    
    partition = "short"
    time      = "12:0:0"
    if not os.path.exists(scratch+geneID+'/'):
        os.mkdir(scratch+geneID+'/')
    output    = scratch+geneID+'/'+ligand_inchikey + "_summary.txt"
    error     = scratch+geneID+'/'+ligand_inchikey + "_error.txt"

    cmd = ' '.join(["sbatch", 
           "-p", partition, 
           "-t", time, 
           "-c", args.ncpus, 
           "-o", output,
           "-e", error, 
           '--exclude="compute-f-17-[09,10-25]"', 
           "--wrap='python " + wd+ "dockpipe_qvina.py {0} {1} {2} {3} {4}'".format(ligand_inchikey, geneID, 
                                                                                   subdomain, str(int(runs)), 
                                                                                   args.ncpus)])          
    print (cmd)
                
"""docking"""
# load the homology models available
runs_dict=pickle.load(open(('../resources/runs_dict_seqid=40_rmsdcutoff=0.1.p'),'rb'))

# main docking 
for index, arow in data.iloc[int(args.start):int(args.end)].iterrows():
    geneID = str(arow['ENTREZ_GENE_ID'])
    ligand_smiles = arow['CANONICAL_SMILES_DESALT']
    ligand_inchikey = arow['INCHIKEY_DESALT']
    
    # first check if target information or homology model is available
    try:
        kinase_domain_name = geneidtokinasename[geneID]
    except KeyError:
        print ('pairs not run because no key for target ', geneID)
        continue
    if not any(kinase_domain_name == fullname[:-3] for fullname in runs_dict.keys()):
        print ('homology models not created for ', kinase_domain_name)
        continue
    
    # run docking if the pair has not yet been docked
    if not os.path.exists(dataset+geneID+'/'+ligand_inchikey+'.zip'):
        lgdtmp = Setup_Fldrs(geneID, ligand_inchikey, ligand_smiles)
        # TODO: assume the first kinase domain is where ligand binds
        # for multi-domain kinase docking into both or need better annotation 
        subdomain = kinase_domain_name + '_D0'
        try:
            runs = Load_RunNumber(lgdtmp, subdomain)
            Submit_Docking(geneID, ligand_inchikey, subdomain, runs)
        except ZeroDivisionError: 
            print ('pairs not run because no homology models after filtering')
            continue
