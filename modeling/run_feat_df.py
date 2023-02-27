#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:13:58 2019

to submit batch jobs to featurize the selected kinase-compound pairs

@author: cl321
"""
import os

datasetdir = '/n/groups/alquraishi/Changchang/Datasets/D13/'
wd = datasetdir + 'M2/I2/'

def Submit_One_Job(start, end):

    partition = 'short'
    time      = '10:0:0'
    output    = '/n/scratch3/users/c/cl321/Datasets/D13/M2/I2/feat/summary/feat_min_pose_'+start+'-'+end+'.txt'
    error     = '/n/scratch3/users/c/cl321/Datasets/D13/M2/I2/feat/error/feat_min_pose_'+start+'-'+end+'.txt'

    cmd = ' '.join(["sbatch", 
           "-p", partition, 
           "-t", time, 
           "-o", output,
           "-e", error, 
           "--mem", '8000',
           '--exclude="compute-f-17-[09,10-25]"', 
           "--wrap='python " + wd+ "feat_dock.py {0} {1}'".format(start, end)])
    print(cmd)

configs = ['cv1', 'cv2', 'cv3', 'cv4', 'cv5']
    
for i in range(0, 283000, 1000):
    for aconf in configs:
        if not os.path.exists(wd+ aconf + '/training_set_' + str(i)+'-'+str(i+1000)+'.hdf'):
            Submit_One_Job(str(i), str(i+1000))
            break
