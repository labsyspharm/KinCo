#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:10:22 2017

This script dispatches jobs to calculate the predicted energy by vina. 

@author: cl321
"""
import pandas as pd
import subprocess
import os
import argparse
import glob
import pickle
import random
import datetime
import sys

datasetdir = '/n/groups/alquraishi/Changchang/Datasets/'
wd = '/n/groups/alquraishi/Changchang/Datasets/vina_energy/'
biochem_df = pd.read_csv("../resources/DTC_cvsplit1127.csv",
                         low_memory = False, dtype={'ENTREZ_GENE_ID':object,
                                               'INCHIKEY_DESALT':object})
biochem_df.drop_duplicates(subset = ['ENTREZ_GENE_ID', 'INCHIKEY_DESALT'], inplace=True)
biochem_df.reset_index(drop=True, inplace=True)

start = int(sys.argv[1])
end = int(sys.argv[2])
batch_size = int(sys.argv[3])

def Submit_One_Job(start, end):

    partition = 'short'
    time      = '5:0:0'
    
    output    = wd + '/tmp/summary/energy_'+start+'-'+end+'.txt'
    error    =  wd + '/tmp/error/energy_'+start+'-'+end+'.txt'

    cmd = ' '.join(["sbatch",
           "-p", partition,
           "-t", time,
           "-o", output,
           "-e", error,
           "--mem", '8000',
           '--exclude="compute-f-17-[09,10-25]"',
           "--wrap='python " + wd+ "vina_rescore.py {0} {1}'".format(start, end)])
    print(cmd)

acc = 0
start_idx = start
for idx, arow in biochem_df.iloc[start:end].iterrows():
    geneid = str(arow['ENTREZ_GENE_ID'])
    lgdinchikey = arow['INCHIKEY_DESALT']
    if not os.path.exists('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+ geneid + '/'+lgdinchikey+'_energy.csv'):
        acc+=1
    if acc == batch_size: 
        acc = 0
        Submit_One_Job(str(start_idx), str(idx+1))
        start_idx = idx+1
if acc!=0:
    Submit_One_Job(str(start_idx), str(idx+1))
