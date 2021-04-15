#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:52:05 2019

After the featurization scripts (feat_dock.py and feat_xtal.py), there are 
many batches of .hdf files. This script will combine them into cvsets for 
training/validation/test. 

Need to tell the script which cv split (1-5) to combine hdfs for

@author: cl321
"""

import sys
import h5py
import glob
import pickle
import numpy as np
import pandas as pd
import os

cvconfig = sys.argv[1]

wd = os.getcwd()+'/%s/'%cvconfig

def read_hdf(filename):
    ids = []
    affinity = []
    coords = []
    features = []
    with h5py.File(filename, 'r') as f:
        for pdb_id in f:
            dataset = f[pdb_id]
            coords.append(dataset[:, :3])
            features.append(dataset[:, 3:])
            affinity.append(dataset.attrs['affinity'])
            ids.append(pdb_id)
    return [ids, affinity, coords, features]

for cvset in ['training', 'validation', 'test']:
    # identify all the batch .hdf files, including those with crystal structures
    filelist = sorted(glob.glob('{0}/{1}_set_*-*.hdf'.format(wd, cvset))) + sorted(glob.glob('{0}/{1}_set_xtal.hdf'.format(wd, cvset)))

    # read each batch .hdf file
    for afile in filelist:
        print ('adding ', afile)
        [ids, affinity, coords, features] = read_hdf(afile)
        
        # add the pairs to the appropriate .hdf file
        with h5py.File(wd+'/%s_set.hdf'%cvset, 'a') as h:

            for index, name in enumerate(ids):
                data = np.hstack([coords[index], features[index]])
                dataset = h.create_dataset(name, data=data, shape=data.shape,
                                       dtype='float32', compression='lzf')
                dataset.attrs['affinity'] = affinity[index]
