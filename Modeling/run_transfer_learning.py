#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Oct 31 16:38:40 2018

@author: cl321
"""
import subprocess
import sys
import os

exp_index = sys.argv[1]
run_index = sys.argv[2]

wd = os.getcwd() 
savedir = wd + '/' + exp_index + '/' + run_index + '/'
computer = wd.split('Research')[0] + '/'
prmtfile = savedir + 'parameters.log'

def fileparser(lines, keyword):
    for aline in lines:
        if keyword in aline:
            result = aline.split()[1:]
            if len(result) == 1:
                return result[0]
            else:
                return result

# read in the parameters in paramters.log
with open(prmtfile, 'r') as pf:
    lines = pf.readlines()
    input_dir = str(fileparser(lines, 'input_dir'))
    lr = str(10**float(fileparser(lines, 'lr')))
    l2 = str(10**float(fileparser(lines, 'l2')))
    random_seed = str(fileparser(lines, 'random_seed'))
    fine_tune = str(fileparser(lines, 'fine_tune')) # False: retrain FC layers; True: freeze FC layers
    meta_path = str(fileparser(lines, 'meta_path'))
    num_epochs = str(fileparser(lines, 'num_epochs'))
    pf.close()

cmd = ' '.join(['python', os.path.dirname(os.path.abspath(__file__)) + '/training_transfer.py',
       '--input_dir', computer + input_dir, '--output_prefix', 
       savedir + '/output', '--log_dir', savedir+'/logdir/', 
       '--l2', l2, '--learning_rate', lr,
       '--num_epochs', num_epochs, '--num_checkpoints', '100', 
       '--fine_tune', fine_tune, '--random_seed', random_seed,
       '--meta_path', computer + meta_path])
print(cmd)
