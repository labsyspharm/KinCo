#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 09:28:19 2017

@author: cl321
"""
import os
import pickle
import numpy as np
import sys
import pandas as pd

parentfldr = os.getcwd()+'/'
seqid = sys.argv[1]
rmsdcutoff = sys.argv[2]
biochem_df = pd.read_csv('/n/groups/alquraishi/Changchang/git_repos/KinCo/resources/the_544_curated_kinome_annotated_v3_short.csv')
geneids = biochem_df['gene_id'].drop_duplicates().astype(str)

geneidtokinasename=pickle.load(open(("./GeneIDtoKinasename.p"),'rb'))

if not os.path.isfile((parentfldr+"/runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".p")):
    runs_dict={}
    not_in_dict = []
    for agene in geneids:
        
        try:
            a_kinase_domain = geneidtokinasename[agene] + '_D0'
        except KeyError:
            print ('gene id not recognized ', agene)
            
        seqiddict={}
        try:
            seqlist = ('/n/groups/alquraishi/Changchang/Ensembler/database_it2/modelrmsds/'
                       + a_kinase_domain + "_sequence-identities.txt")
            with open(seqlist,"r") as f:
                for line in f:
                    linesplit=line.split()
                    seqiddict[linesplit[0]]=linesplit[1]
            with open(('/n/groups/alquraishi/Changchang/Ensembler/database_it2/modelrmsds/'
                + a_kinase_domain + "_" + rmsdcutoff + "_unique-models.txt"),"r") as g:
                uniquemodels=g.read().splitlines()
                if 'KPCB_RAT_D0_3PFQ_A' in uniquemodels:
                    uniquemodels.remove('KPCB_RAT_D0_3PFQ_A')
                filtered_models = []
                for amodel in uniquemodels:
                    try:
                        if float(seqiddict[amodel]) >= int(seqid):
                        		filtered_models.append(amodel)
                    except KeyError:
                        print (amodel, a_kinase_domain)
            runs_dict[a_kinase_domain] = filtered_models
            if len(filtered_models) == 0:
                print ("Warning: no models produced ", a_kinase_domain)
                not_in_dict.append(agene)
        except IOError: 
            continue
    with open((parentfldr+"/runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".p"),"wb") as handle:
                pickle.dump(runs_dict,handle)
    np.savetxt(parentfldr+"/not_in_runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".txt", not_in_dict, fmt='%s')
else:
    print ('loaidng existing runs_dict')
    
    runs_dict=pickle.load(open((parentfldr+"/runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".p"),'rb'))
    not_in_dict = list(np.loadtxt(parentfldr+"/not_in_runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".txt", dtype=str))
    for agene in geneids:

        try:
            a_kinase_domain = geneidtokinasename[agene] + '_D0'
        except KeyError:
            print ('gene id not recognized ', agene)
            
        seqiddict={}
        try:
            seqlist = ('/n/groups/alquraishi/Changchang/Ensembler/database_it2/modelrmsds/'
                       + a_kinase_domain + "_sequence-identities.txt")
            with open(seqlist,"r") as f:
                for line in f:
                    linesplit=line.split()
                    seqiddict[linesplit[0]]=linesplit[1]
            with open(('/n/groups/alquraishi/Changchang/Ensembler/database_it2/modelrmsds/'
                + a_kinase_domain + "_" + rmsdcutoff + "_unique-models.txt"),"r") as g:
                uniquemodels=g.read().splitlines()
                if 'KPCB_RAT_D0_3PFQ_A' in uniquemodels:
                    uniquemodels.remove('KPCB_RAT_D0_3PFQ_A')
                filtered_models = []
                for amodel in uniquemodels:
                    try:
                        if float(seqiddict[amodel]) >= int(seqid):
                            filtered_models.append(amodel)
                    except KeyError:
                        print (amodel, a_kinase_domain)
            runs_dict[a_kinase_domain] = filtered_models
            if len(filtered_models) == 0:
                print ("Warning: no models produced ", a_kinase_domain)
                not_in_dict.append(agene)
        except IOError: 
            continue
    with open((parentfldr+"/runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".p"),"wb") as handle:
                pickle.dump(runs_dict,handle)
    np.savetxt(parentfldr+"/not_in_runs_dict_seqid="+seqid+"_rmsdcutoff="+rmsdcutoff+".txt", not_in_dict, fmt='%s')
