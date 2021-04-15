""" 

script adapted from https://gitlab.com/cheminfIBB/pafnucy/training.py 
publication: Marta M Stepniewska-Dziubinska, Piotr Zielenkiewicz, Pawel Siedlecki, 
Development and evaluation of a deep learning model for protein–ligand binding 
affinity prediction, Bioinformatics, Volume 34, Issue 21, 01 November 2018, 
Pages 3666–3674, https://doi.org/10.1093/bioinformatics/bty374

The script trains 3DCNN using all PDBBIND2018 data. MSE is used in the loss function
to be run with run_pretraining.py

"""

import numpy as np

import pandas as pd
from math import sqrt, ceil

import h5py

from sklearn.utils import shuffle
import tensorflow as tf

from tfbio.data import Featurizer, make_grid, rotate
import tfbio.net_noReLU as net_noReLU

import os.path

import matplotlib as mpl
mpl.use('agg')


import time
timestamp = time.strftime('%Y-%m-%dT%H:%M:%S')


datasets = ['training', 'validation', 'test']


def input_dir(path):
    """Check if input directory exists and contains all needed files"""
    global datasets

    path = os.path.abspath(path)
    if not os.path.isdir(path):
        raise IOError('Incorrect input_dir specified: no such directory')
    for dataset_name in datasets:
        dataset_path = os.path.join(path, '%s_set.hdf' % dataset_name)
        if not os.path.exists(dataset_path):
            raise IOError('Incorrect input_dir specified:'
                          ' %s set file not found' % dataset_path)
    return path

import argparse
parser = argparse.ArgumentParser(
    description='Train 3D colnvolutional neural network on affinity data',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

io_group = parser.add_argument_group('I/O')
io_group.add_argument('--input_dir', '-i', required=True, type=input_dir,
                      help='directory with training, validation and test sets')
io_group.add_argument('--log_dir', '-l', default='./logdir/',
                      help='directory to store tensorboard summaries')
io_group.add_argument('--output_prefix', '-o', default='./output',
                      help='prefix for checkpoints, predictions and plots')
io_group.add_argument('--grid_spacing', '-g', default=1.0, type=float,
                      help='distance between grid points')
io_group.add_argument('--max_dist', '-d', default=10.0, type=float,
                      help='max distance from complex center')

arc_group = parser.add_argument_group('Netwrok architecture')
arc_group.add_argument('--conv_patch', default=5, type=int,
                       help='patch size for convolutional layers')
arc_group.add_argument('--pool_patch', default=2, type=int,
                       help='patch size for pooling layers')
arc_group.add_argument('--conv_channels', metavar='C', default=[64, 128, 256],
                       type=int, nargs='+',
                       help='number of fileters in convolutional layers')
arc_group.add_argument('--dense_sizes', metavar='D', default=[1000, 500, 200],
                       type=int, nargs='+',
                       help='number of neurons in dense layers')

reg_group = parser.add_argument_group('Regularization')
reg_group.add_argument('--keep_prob', dest='kp', default=0.5, type=float,
                       help='keep probability for dropout')
reg_group.add_argument('--l2', dest='lmbda', default=0.001, type=float,
                       help='lambda for weight decay')
reg_group.add_argument('--rotations', metavar='R', default=list(range(24)),
                       type=int, nargs='+',
                       help='rotations to perform')

tr_group = parser.add_argument_group('Training')
tr_group.add_argument('--learning_rate', default=1e-5, type=float,
                      help='learning rate')
tr_group.add_argument('--batch_size', default=20, type=int,
                      help='batch size')
tr_group.add_argument('--num_epochs', default=20, type=int,
                      help='number of epochs')
tr_group.add_argument('--num_checkpoints', dest='to_keep', default=10, type=int,
                      help='number of checkpoints to keep')
tr_group.add_argument('--random_seed', default=123, type=int,
                      help='random seed to control stochasticity')
args = parser.parse_args()

prefix = os.path.abspath(args.output_prefix) + '-' + timestamp
logdir = os.path.join(os.path.abspath(args.log_dir), os.path.split(prefix)[1])
np.random.seed(args.random_seed)
print(args.output_prefix[:-6])
print(prefix)

# load and process data
featurizer = Featurizer()

print('\n---- FEATURES ----\n')
print('atomic properties:', featurizer.FEATURE_NAMES)

columns = {name: i for i, name in enumerate(featurizer.FEATURE_NAMES)}

ids = {}
affinity = {}
coords = {}
features = {}

for dictionary in [ids, affinity, coords, features]:
    for dataset_name in datasets:
        dictionary[dataset_name] = []

for dataset_name in datasets:
    dataset_path = os.path.join(args.input_dir, '%s_set.hdf' % dataset_name)
    with h5py.File(dataset_path, 'r') as f:
        for pdb_id in f:
            dataset = f[pdb_id]

            coords[dataset_name].append(dataset[:, :3])
            features[dataset_name].append(dataset[:, 3:])
            affinity[dataset_name].append(dataset.attrs['affinity'])
            ids[dataset_name].append(pdb_id)

    ids[dataset_name] = np.array(ids[dataset_name])
    affinity[dataset_name] = np.reshape(affinity[dataset_name], (-1, 1))


# normalize charges
charges = []
for feature_data in features['training']:
    charges.append(feature_data[..., columns['partialcharge']])

charges = np.concatenate([c.flatten() for c in charges])

m = charges.mean()
std = charges.std()
print('charges: mean=%s, sd=%s' % (m, std))
print('use sd as scaling factor')
if not os.path.exists(args.input_dir + '/scaling_factor.txt'):
    np.savetxt(args.input_dir + '/scaling_factor.txt', [np.round(std, 3)])

def get_batch(dataset_name, indices, rotation=0):
    global coords, features, std
    if dataset_name == 'training':
        std_round = std
    else:
        std_round = np.round(std, 3)

    x = []
    for i, idx in enumerate(indices):
        coords_idx = rotate(coords[dataset_name][idx], rotation)
        features_idx = features[dataset_name][idx]
        x.append(make_grid(coords_idx, features_idx,
                 grid_resolution=args.grid_spacing,
                 max_dist=args.max_dist))
    x = np.vstack(x)
    x[..., columns['partialcharge']] /= std_round
    return x

def get_predictions(val_or_test):

    rot_pred = {}
    x_v = range(ds_sizes[val_or_test])
    
    for arot in args.rotations:
        arot_pred = []
        for b, (bi, bj) in enumerate(batches(val_or_test)):
            
            arot_pred.append(session.run(y, feed_dict={x: get_batch(val_or_test, x_v[bi:bj], arot),
                       keep_prob: 1.0}))
        rot_pred[str(arot)] = np.vstack(arot_pred).flatten()

    rot_pred_df = pd.DataFrame(rot_pred)
    rot_pred_df['mean'] = rot_pred_df.mean(axis=1)
    rot_pred_df['expt'] = list(affinity[val_or_test].flatten())
    rot_pred_df['name'] = list(ids[val_or_test].flatten())
    rot_pred_df['hdf_id'] = rot_pred_df['name'].str.split('_').str[0]
    merge = rot_pred_df.merge(ref_cvset.drop_duplicates(subset=['hdf_id']), on=['hdf_id'], how  = 'left')
    merge['dataset'] = merge['name'].str.split('_', expand=True)[1]
    merge['RESULT_VALUE_min'] = merge['expt']
    merge.to_csv('{0}/{1}_{2}_prediction.csv'.format(args.output_prefix[:-6], val_or_test, epoch),
                 index=False)
        
    return merge


print('\n---- DATA ----\n')

tmp = get_batch('training', range(min(50, len(features['training']))))

assert ((tmp[:, :, :, :, columns['molcode']] == 0.0).any()
        and (tmp[:, :, :, :, columns['molcode']] == 1.0).any()
        and (tmp[:, :, :, :, columns['molcode']] == -1.0).any()).all()

idx1 = [[i[0]] for i in np.where(tmp[:, :, :, :, columns['molcode']] == 1.0)]
idx2 = [[i[0]] for i in np.where(tmp[:, :, :, :, columns['molcode']] == -1.0)]

print('\nexamples:')
for mtype, mol in [['ligand', tmp[idx1]], ['protein', tmp[idx2]]]:
    print(' ', mtype)
    for name, num in columns.items():
        print('  ', name, mol[0, num])
    print('')


# Best error we can get without any training (MSE from training set mean):
t_baseline = ((affinity['training'] - affinity['training'].mean()) ** 2.0).mean()
print('baseline mse: training=%s' % t_baseline)


# NET PARAMS

ds_sizes = {dataset: len(affinity[dataset]) for dataset in datasets}
_, isize, *_, in_chnls = get_batch('training', [0]).shape
osize = 1

for set_name, set_size in ds_sizes.items():
    print('%s %s samples' % (set_size, set_name))

num_batches = {dataset: ceil(size / args.batch_size)
               for dataset, size in ds_sizes.items()}

print('\n---- MODEL ----\n')
print((isize - 1) * args.grid_spacing, 'A box')
print(in_chnls, 'features')
print('')
print('convolutional layers: %s channels, %sA patch + max pooling with %sA patch'
      % (', '.join((str(i) for i in args.conv_channels)), args.conv_patch,
         args.pool_patch))
print('fully connected layers:', ', '.join((str(i) for i in args.dense_sizes)),
      'neurons')
print('regularization: dropout (keep %s) and L2 (lambda %s)'
      % (args.kp, args.lmbda))
print('')
print('learning rate', args.learning_rate)
print(num_batches['training'], 'batches,', args.batch_size, 'examples each')
print('')
print(args.num_epochs, 'epochs, best', args.to_keep, 'saved')

config=tf.ConfigProto(inter_op_parallelism_threads=2, 
                      intra_op_parallelism_threads=2)
config.gpu_options.allow_growth = True
#config.gpu_options.per_process_gpu_memory_fraction = 0.5

# set up the reference experimental table for prediction
ref_cvset = pd.read_csv('/data/cl321/Datasets/PDBBind2018/M2/I9/PDBBIND_cvsplit1127.csv',
        dtype={'ENTREZ_GENE_ID':object})
cvconfig = 'cv' + str(args.input_dir.split('cv')[-1][:1])
ref_cvset['cthresh'] = ref_cvset[cvconfig+'|cthresh']
ref_cvset['tthresh'] = ref_cvset[cvconfig+'|tthresh']
ref_cvset['cvset'] = ref_cvset[cvconfig+'|cvset'] 

# import network
graph = net_noReLU.make_SB_network(isize=isize, in_chnls=in_chnls, osize=osize,
                                  conv_patch=args.conv_patch,
                                  pool_patch=args.pool_patch,
                                  conv_channels=args.conv_channels,
                                  dense_sizes=args.dense_sizes,
                                  lmbda=args.lmbda,
                                  learning_rate=args.learning_rate,
                                  seed = args.random_seed)


train_writer = tf.summary.FileWriter(os.path.join(logdir, 'training_set'),
                                     graph, flush_secs=1)

net_summaries, training_summaries = net_noReLU.make_summaries_SB(graph)

x = graph.get_tensor_by_name('input/structure:0')
y = graph.get_tensor_by_name('output/prediction:0')
t = graph.get_tensor_by_name('input/affinity:0')
keep_prob = graph.get_tensor_by_name('fully_connected/keep_prob:0')
train = graph.get_tensor_by_name('training/train:0')
mse = graph.get_tensor_by_name('training/mse:0')
feature_importance = graph.get_tensor_by_name('net_properties/feature_importance:0')
global_step = graph.get_tensor_by_name('training/global_step:0')

convs = '_'.join((str(i) for i in args.conv_channels))
fcs = '_'.join((str(i) for i in args.dense_sizes))

with graph.as_default():
    saver = tf.train.Saver(max_to_keep=args.to_keep)


def batches(set_name):
    """Batch generator, yields slice indices"""
    global num_batches, args, ds_sizes
    for b in range(num_batches[set_name]):
        bi = b * args.batch_size
        bj = (b + 1) * args.batch_size
        if b == num_batches[set_name] - 1:
            bj = ds_sizes[set_name]
        yield bi, bj

err = float('inf')

train_sample = min(args.batch_size, len(features['training']))

# model training
print('\n---- TRAINING ----\n')
with tf.Session(graph=graph, config = config) as session:
    session.run(tf.global_variables_initializer())
    print (args.output_prefix[:-6])
    latest_ckpt = tf.train.latest_checkpoint(args.output_prefix[:-6])
    if latest_ckpt != None:
        print ('loaded check point ', latest_ckpt)
        saver.restore(session, latest_ckpt)
        current_epoch = int(int(latest_ckpt.split('-')[-1])/(24*int(num_batches['training'])))
    else:
        current_epoch = 0
    print('current epoch', current_epoch)
        
    summary_imp = tf.Summary()
    feature_imp = session.run(feature_importance)
    image = net_noReLU.feature_importance_plot(feature_imp)
    summary_imp.value.add(tag='feature_importance_%s' % 0, image=image)
    train_writer.add_summary(summary_imp, 0)

    stats_net = session.run(
        net_summaries,
        feed_dict={x: get_batch('training', range(train_sample)),
                   t: affinity['training'][:train_sample],
                   keep_prob: 1.0}
    )

    train_writer.add_summary(stats_net, 0)

    for epoch in range(current_epoch, args.num_epochs):
        for rotation in args.rotations:
            print('rotation', rotation)
            # TRAIN #
            x_t, y_t = shuffle(range(ds_sizes['training']), affinity['training'])

            for bi, bj in batches('training'):
                session.run(train, feed_dict={x: get_batch('training',
                                                           x_t[bi:bj],
                                                           rotation),
                                              t: y_t[bi:bj], keep_prob: args.kp})

            # SAVE STATS - per rotation #
            stats_t, stats_net = session.run(
                [training_summaries, net_summaries],
                feed_dict={x: get_batch('training', x_t[:train_sample]),
                           t: y_t[:train_sample],
                           keep_prob: 1.0}
            )

            train_writer.add_summary(stats_t, global_step.eval())
            train_writer.add_summary(stats_net, global_step.eval())

        # SAVE STATS - per epoch #
        # training set error
        pred_t = np.zeros((ds_sizes['training'], 1))
        mse_t = np.zeros(num_batches['training'])

        for b, (bi, bj) in enumerate(batches('training')):
            weight = (bj - bi) / ds_sizes['training']

            pred_t[bi:bj], mse_t[b] = session.run(
                [y, mse],
                feed_dict={x: get_batch('training', x_t[bi:bj]),
                           t: y_t[bi:bj],
                           keep_prob: 1.0}
            )

            mse_t[b] *= weight

        mse_t = mse_t.sum()

        summary_mse = tf.Summary()
        summary_mse.value.add(tag='mse_all', simple_value=mse_t)
        train_writer.add_summary(summary_mse, global_step.eval())

        # predictions distribution
        summary_pred = tf.Summary()
        summary_pred.value.add(tag='predictions_all',
                               histo=net_noReLU.custom_summary_histogram(pred_t))
        train_writer.add_summary(summary_pred, global_step.eval())

        # SAVE MODEL #
        print('epoch: %s train error: %s' % (epoch, mse_t))
        
        checkpoint = saver.save(session, prefix, global_step=global_step)
        
        # validation set error
        get_predictions('validation')
        # test set error
        get_predictions('test')
        
train_writer.close()
