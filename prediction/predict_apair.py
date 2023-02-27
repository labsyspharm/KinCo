#!/usr/bin/env python3
""" 
This script predicts the affinities of kinase-compound pairs as specified in 
input. Due to storage space constraint, I will pull the hdf from the cloud 
and make predictions locally (instead of generating the hdf locally or have
pretransfer them ahead of time). 

script adapted from https://gitlab.com/cheminfIBB/pafnucy/predict.py 
publication: Marta M Stepniewska-Dziubinska, Piotr Zielenkiewicz, Pawel Siedlecki, 
Development and evaluation of a deep learning model for protein–ligand binding 
affinity prediction, Bioinformatics, Volume 34, Issue 21, 01 November 2018, 
Pages 3666–3674, https://doi.org/10.1093/bioinformatics/bty374
To speed up the process, the grid generating step is re-implemented in TensorFlow

"""
import time
total_begin = time.time()
import random

import numpy as np
import pandas as pd
import h5py

import tensorflow as tf
from tensorflow.python.client import timeline

from tfbio.data import Featurizer, rotation_matrix
from tfbio.net_noReLU_v1 import convolve3D, feedforward

from math import ceil, pi
from itertools import combinations

import os
import glob

import subprocess as sp
import shutil

def fileparser(lines, keyword):
    for aline in lines:
        if keyword in aline:
            result = aline.split()[1:]
            if len(result) == 1:
                return result[0]
            else:
                return result
            
def input_file(path):
    """Check if input file exists."""

    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise IOError('File %s does not exist.' % path)
    return path


def network_prefix(paths):
    """Check if all file required to restore the network exists."""
    from glob import glob
    
    for path in paths.split(','):
        dir_path, file_name = os.path.split(path)
        path = os.path.join(os.path.abspath(dir_path), file_name)
    
        for extension in ['index', 'meta', 'data*']:
            file_name = '%s.%s' % (path[:-5], extension)
    
            # use glob instead of os because we need to expand the wildcard
            if len(glob(file_name)) == 0:
                raise IOError('File %s does not exist.' % file_name)

    return paths


def batch_size(value):
    """Check if batch size is a non-negative integer"""

    value = int(value)
    if value < 0:
        raise ValueError('Batch size must be positive, %s given' % value)
    return value


def output_file(path):
    """Check if output file can be created."""

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    if not os.access(dirname, os.W_OK):
        raise IOError('File %s cannot be created (check your permissions).'
                      % path)
    return path


def string_bool(s):
    s = s.lower()
    if s in ['true', 't', '1', 'yes', 'y']:
        return True
    elif s in ['false', 'f', '0', 'no', 'n']:
        return False
    else:
        raise IOError('%s cannot be interpreted as a boolean' % s)
        
def read_hdf(input_hdf):
    """
    to read in a hdf (all poses of a kinase-compoun dpair), load into memory
    """
    featurizer = Featurizer()
    
    charge_column = featurizer.FEATURE_NAMES.index('partialcharge')
    
    coords = []
    features = []
    names = []
    
    with h5py.File(input_hdf, 'r') as f:
        for name in f:
            names.append(name)
            dataset = f[name]
            coords.append(dataset[:, :3])
            features.append(dataset[:, 3:])
    
    coords = np.asarray(coords)
    features = np.asarray(features)
    names = np.asarray(names)
    
    return [coords, features, names, charge_column]

def tf_rotate(coords, rotation):
    """Rotate coordinates by a given rotation
    """
    global ROTATIONS
    
    ROTATIONS_tf = tf.constant(np.asarray(ROTATIONS))
    rot_mtx_tf = tf.tile(tf.expand_dims(ROTATIONS_tf[rotation],0), [tf.shape(ROTATIONS_tf[rotation])[0], 1, 1])
    return tf.matmul(coords, tf.cast(rot_mtx_tf, tf.float32))

def tf_put_box(grid_coords, box_size = 21):
    """select on the coordinates that fit within the box"""
    #coordinates less than 10 A is less than box size
    mask1 = tf.less(grid_coords, box_size)
    #coordinates greater than -10 A is greater than 0
    mask2 = tf.greater_equal(grid_coords,0)
    in_box_bool = tf.reduce_all(mask1&mask2, axis = 2)
    
    return in_box_bool

def tf_make_grid(coords_rotated, features, sel_atoms, labels):
    num_features = 19
    box_size = 21
    new_coords = tf.cast(tf.boolean_mask(coords_rotated,sel_atoms), tf.int32)
    new_features = tf.boolean_mask(features, sel_atoms)
    filled_grid = tf.scatter_nd(new_coords, new_features, 
                                shape=[box_size, box_size, box_size, num_features])
    return [filled_grid, labels]


def make_pred_network(graph, network_input, isize=20, in_chnls=19, osize=1,
                    conv_patch=5, pool_patch=2, conv_channels=[64, 128, 256],
                    dense_sizes=[1000, 500, 200],
                    lmbda=0.001, learning_rate=1e-5,
                    seed=123):
    
    with graph.as_default(): 
        with tf.variable_scope('pafnucy/convolution'):
            h_convs = convolve3D(network_input, conv_channels,
                                 conv_patch=conv_patch,
                                 pool_patch=pool_patch)
            
        hfsize = isize
        for _ in range(len(conv_channels)):
            hfsize = ceil(hfsize / pool_patch)
        hfsize = conv_channels[-1] * hfsize**3

        with tf.variable_scope('transfer/fully_connected'):
            h_flat = tf.reshape(h_convs, shape=(-1, hfsize), name='h_flat')

            prob1 = tf.constant(1.0, name='keep_prob_default')
            keep_prob = tf.placeholder_with_default(prob1, shape=(),
                                                    name='keep_prob')

            h_fcl = feedforward(h_flat, dense_sizes, keep_prob=keep_prob)

        with tf.variable_scope('transfer/output'):
            w = tf.get_variable('w', shape=(dense_sizes[-1], osize),
                                initializer=tf.truncated_normal_initializer(
                                    stddev=(1 / (dense_sizes[-1]**0.5))))
            b = tf.get_variable('b', shape=(osize,), dtype=tf.float32,
                                initializer=tf.constant_initializer(1))
            y = tf.add(tf.matmul(h_fcl, w), b, name='prediction')

    return graph

# Create matrices for all possible 90 degree rotations of a box
ROTATIONS = [rotation_matrix([1, 1, 1], 0)]
# about X, Y and Z - 9 rotations
for a1 in range(3):
    for t in range(1, 4):
        axis = np.zeros(3)
        axis[a1] = 1
        theta = t * pi / 2.0
        ROTATIONS.append(rotation_matrix(axis, theta))
# about each face diagonal - 6 rotations
for (a1, a2) in combinations(range(3), 2):
    axis = np.zeros(3)
    axis[[a1, a2]] = 1.0
    theta = pi
    ROTATIONS.append(rotation_matrix(axis, theta))
    axis[a2] = -1.0
    ROTATIONS.append(rotation_matrix(axis, theta))
# about each space diagonal - 8 rotations
for t in [1, 2]:
    theta = t * 2 * pi / 3
    axis = np.ones(3)
    ROTATIONS.append(rotation_matrix(axis, theta))
    for a1 in range(3):
        axis = np.ones(3)
        axis[a1] = -1
        ROTATIONS.append(rotation_matrix(axis, theta))
        
def predict_affinity(coords, features, names, charge_column, network, 
                     grid_spacing, 
                     rotations, 
                     max_dist,
                     batch, 
                     verbose = True):
    
    """
    predict the binding affinity for a kinase-compound pair (one hdf file,
    multiple poses) based on a model with indicated number of rotations
    
    input: 
        coords: the coordinates of all poses (from hdf)
        features: the features of all poses (from hdf)
        names: the poseID and homology of the pose (from hdf)
        charge_column: the column with the partial charge
        network: path of the model to use
        grid_spacing: resolution of the grid
        rotations: rotation matrix. How many times to rotate the complex
        max_dist: distance from the center of the ligand
        batch: number of poses to process at the same time
    output: 
        a table with the predicted affinity for all poses of the kinase-compound pair 
        (all complexes in the hdf)
    """
    global ROTATIONS
    
    if args.verbose:
        print('loaded %s complexes\n' % len(coords))

    config = tf.ConfigProto()
    config.intra_op_parallelism_threads = 4
    config.inter_op_parallelism_threads = 4
    
    tf.reset_default_graph()
    
    graph = tf.Graph()

    # modify the graph for prediction by adding grid featurization 
    # and rotation for speed up
    with graph.as_default():
        features_placeholder = tf.placeholder(features.dtype, features.shape)
        coords_placeholder = tf.placeholder(coords.dtype, coords.shape)
        labels_placeholder = tf.placeholder(names.dtype, names.shape)
        
        # featurize the coordinates into grid box (including rotation) on the graph
        rotation_placeholder = tf.placeholder(tf.int32)
        ROTATIONS_tf = tf.constant(np.asarray(ROTATIONS))
        rot_mtx = tf.tile(tf.expand_dims(ROTATIONS_tf[rotation_placeholder],0), 
                          [tf.shape(coords_placeholder)[0], 1, 1])
        
        rot_coords = tf.matmul(tf.cast(coords_placeholder, tf.float32),
                               tf.cast(rot_mtx, tf.float32))
        grid_coords = (rot_coords + max_dist) / grid_spacing
        grid_coords = tf.cast(tf.round(grid_coords), tf.int32)

        in_box_bool = tf_put_box(grid_coords)

        # set up data queue
        dataset = tf.data.Dataset.from_tensor_slices((grid_coords, 
                                              features_placeholder, 
                                              in_box_bool,
                                              labels_placeholder))
    
        dataset = dataset.map(lambda c, f, b, l: tf_make_grid(c, f, b, l))

        dataset = dataset.batch(args.batch)
        dataset = dataset.prefetch(buffer_size=args.batch)
        iterator = dataset.make_initializable_iterator()
        [grid, label_input] = iterator.get_next()
    
        # load the model network
        graph = make_pred_network(graph, grid)
    
    # get the prediction fromt he graph
    y = graph.get_tensor_by_name('transfer/output/prediction:0')
    # get the parameters from the graph to retrieve trained values
    with graph.as_default():
        var_extract = tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope='')
        
    var_dict = {}
    for var in var_extract:
        var_dict[var.name] = var

    pairnames = []
    predictions = []
    results = {}
    
    with tf.Session(graph=graph, config=config) as session:
        
        # predict affinity for each pose and each rotation
        for index, i in enumerate(args.rotations):
            predictions = []
            session.run(iterator.initializer, feed_dict={features_placeholder: features,
                                                        coords_placeholder: coords,
                                                        labels_placeholder: names,
                                                        rotation_placeholder: i})
            
            pafnucy_saver = tf.train.Saver(var_list = var_extract)
            pafnucy_saver.restore(session, network[:-5])
                
            while True:
                try:
                    if index == 0:
                        [predictions_batch, names_batch] = session.run([y, label_input])
                        predictions+=list(np.reshape(predictions_batch, -1))
                        pairnames+=list(np.reshape(names_batch, -1))
                    else:
                        predictions_batch = session.run(y)
                        predictions+=list(np.reshape(predictions_batch, -1))
                except tf.errors.OutOfRangeError:
                    break
      
            results[str(i)] = predictions

    # process ID and predictions into table
    names_split = [i.split('|') for i in names]  
    modelnames = [i[0] for i in names_split]
    poseIDs = [i[1] for i in names_split]
    
    results['modelname'] = modelnames 
    results['poseID'] = poseIDs
    results = pd.DataFrame(results)
    
    if len(args.rotations) != 1:
        results['Kdmodel'] = results.mean(axis=1)
    else:
        results.rename(columns={'0':'Kdmodel'}, inplace=True)
            
    return results

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(
        description='Predict affinity with the network',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''This script reads the structures of complexes from HDF file and
        predicts binding affinity for each comples. 
        '''
    )
    
    parser.add_argument('--inputfile', '-i', required=True, type=str,
                        help='file with target gene IDs and ligand inchikeys')
    parser.add_argument('--modelfile', '-m', type=str,
                        help='file with locations of necessary files')
    parser.add_argument('--grid_spacing', '-g', default=1.0, type=float,
                        help='distance between grid points used during training')
    parser.add_argument('--rotations', metavar='R', default=[0],
                           type=int, nargs='+',
                           help='rotations to perform')
    parser.add_argument('--max_dist', '-d', default=10.0, type=float,
                        help='max distance from complex center used during training')
    parser.add_argument('--batch', '-b', type=batch_size,
                        default=20,
                        help='batch size. If set to 0, predict for all complexes at once.')
    parser.add_argument('--verbose', '-v', type=string_bool,
                        default=True,
                        help='whether to print messages')
    
    args = parser.parse_args()
    
    # read in info about which models to use for prediction and the scaling factors for the partial charge
    with open(args.modelfile, 'r') as pf:
        lines = pf.readlines()
        itdir = str(fileparser(lines, 'itdir'))
        iteration = str(fileparser(lines, 'iteration'))
        modelid = str(fileparser(lines, 'modelid')).split(',')
        randomseed = int(fileparser(lines, 'randomseed'))
        random.seed(randomseed)
        scaling_factors = fileparser(lines, 'scaling_factors').split(',')
        modelpaths = str(fileparser(lines, 'modelpaths')).split(',')
        
    datasetdir =  "/data/cl321/Datasets/D13/M2/%s/"%itdir  
    modelparentdir =  "/data/cl321/Datasets/D13/M2/%s/models/"%itdir    
    datasetdir_O2 =  "/n/groups/alquraishi/Changchang/Datasets/D13/M2/%s/"%itdir 
    datasetdir_scratch = '/n/scratch3/users/c/cl321/Datasets/D13/M2/%s/'%itdir
    wd = os.getcwd()+'/'

    myhost = os.uname()[1]
    print('script is running at ', myhost)
    [machine, gpu] = args.inputfile.split('/')[-3:-1]
    
    # read in which kinase-compound pairs to predict 
    pair_list = np.loadtxt(args.inputfile, dtype=str)
    
    for idx, apair in enumerate(pair_list):
        [agene, algd] = apair.split('|')
        print(idx, agene, algd)
        
        inputpathtgt=datasetdir+'/CNN_features/' + agene
        inputpath = inputpathtgt + '/' + algd + '.hdf'
        # locate the HDF feature files for each pair
        inputpath_O2 = '{0}/CNN_features/{1}/{2}/{3}/{4}.hdf'.format(datasetdir_scratch,
                                                        machine, gpu, agene, algd)
        
        # predict affinity if energy file doesn't already exist
        energy_pass = len(glob.glob(datasetdir + '/energy/*/'+agene+'/'+algd+'_tf.csv'))
        if not os.path.exists(inputpath):
            if energy_pass >= len(modelid):
                continue
            else:
                if not os.path.exists(inputpathtgt):
                    os.mkdir(inputpathtgt)
                cmd = ['rsync', '--remove-source-files','-avz',
                       'cl321@transfer.rc.hms.harvard.edu:'+inputpath_O2, inputpath]
                sp.call(cmd)
        
        # read and load features. raise error if something wrong with hdf file
        try:
            [coords, features_raw, names, charge_column] = read_hdf(inputpath)
        except OSError:
            print('no hdf', inputpath)
            with open('error_log.txt', 'a+') as f:
                f.write('OS '+inputpath+'\n')
            if os.path.exists(inputpath):
                os.remove(inputpath) #due to OSError can't read data. will remove and regenerate
            continue
        except RuntimeError:
            with open('error_log.txt', 'a+') as f:
                f.write('RT '+inputpath+'\n')
            os.remove(inputpath)
            continue
        except KeyError:
            os.remove(inputpath)
            with open('error_log.txt', 'a+') as f:
                f.write('KY '+inputpath+'\n')
            print('key error', inputpath)
        
        # predict affinity for each model 
        for midx in range(len(modelid)):
            aNID = modelid[midx]
            anet = modelparentdir + modelpaths[midx]
            charge_scaler = float(scaling_factors[midx])
            features = features_raw.copy()
            features[:, :, 12] /= charge_scaler 
            
            sfsenergy = datasetdir + '/energy/It'+iteration+"_"+ aNID + '/'

            if not os.path.exists(sfsenergy+agene):
                os.mkdir(sfsenergy+agene)
            energyfilepath = sfsenergy+agene+'/'+algd+'_tf.csv'
        

            if not os.path.exists(energyfilepath):
                results =  predict_affinity(coords, features, names, 
                                            charge_column, 
                                            anet,
                                            grid_spacing = args.grid_spacing,
                                            rotations = args.rotations,
                                            max_dist = args.max_dist,
                                            batch = args.batch)
                results.to_csv(energyfilepath, index=False)
        
        for aNID in modelid:
            sfsenergy = datasetdir + '/energy/It'+iteration+"_"+ aNID + '/'
            energyfilepath = sfsenergy+agene+'/'+algd+'_tf.csv'
            delfile = True
            if not os.path.exists(energyfilepath):
                delfile = False
        if delfile:
            os.remove(inputpath)   
        
        

