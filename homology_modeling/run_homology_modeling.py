#!/n/app/python/2.7.12-ucs4/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:33:58 2016

This script dispatches individual homology modeling jobs for each kinase in the
targets/target.fa file. Runs homolgy_modeling.py which runs the
Ensembler package https://ensembler.readthedocs.io/en/latest/

@author: cl321
"""

targets=[]

with open('/n/groups/alquraishi/Changchang/Ensembler/database/targets/targets.fa') as f:
    lines=f.readlines()
    for line in lines:
        if line.startswith('>'):
            name=line.split('>')[1]
            targets.append(name)
            
for i in targets:
    
    ID = i.strip('\n')

    print(' '.join(['sbatch', '-p', 'long', '-t', '240:00:00', 
                    '-o', '/n/scratch3/users/c/cl321/Ensembler/database/summary/summary_'+str(ID), 
                    '-e', '/n/scratch3/users/c/cl321/Ensembler/database/error/error_'+str(ID), 
                    'homology_modeling.py', ID]))
