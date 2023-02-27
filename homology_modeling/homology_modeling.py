#!/n/app/python/2.7.12-ucs4/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:42:52 2016

This script runs the homology modeling for a target kinase using all the template
structures indicated under templates/

It is a wrapper around the Ensembler package https://ensembler.readthedocs.io/en/latest/
and has 3 steps
1) sequence alignment between the target kinase and the template kinases
2) homology modeling based on the sequence alignment
3) cluster the homology models by structure similarity to reduce the number of 
redundant structures

At the end of the pipeline, the ID of the unique homology models and the sequence
identity of the templates to the target kinase will be copied to modelrmsd 
folder to be used as info for quality check and selection of receptors for docking

This script takes the ID of an individual kinase (uniprot name) to start the pipeline

@author: cl321
"""
import argparse
import subprocess
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('ID', help='kinase to process')
args = parser.parse_args()

cutoff = '0.1'

subprocess.call(['ensembler', 'align', '--targets', args.ID])
subprocess.call(['ensembler', 'build_models', '--targets', args.ID])
subprocess.call(['ensembler', 'cluster', '--cutoff', cutoff, '--targets', args.ID])

source1=('/n/groups/alquraishi/Changchang/Ensembler/database/models/'+args.ID+
        '/unique-models.txt')
dest1=('/n/groups/alquraishi/Changchang/Ensembler/database/modelrmsds/'+
        args.ID+'_'+cutoff+'_unique-models.txt')
shutil.copyfile(source1, dest1)
source2=('/n/groups/alquraishi/Changchang/Ensembler/database/models/'+args.ID+
        '/sequence-identities.txt')
dest2=('/n/groups/alquraishi/Changchang/Ensembler/database/modelrmsds/'+
        args.ID+'_sequence-identities.txt')
shutil.copyfile(source2, dest2)
