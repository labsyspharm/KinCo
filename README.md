# KinCo

KinCo is a dataset containing in silico structures of kinase-compound pairs. 
 <br><br>
 <br><br>
 <br><br>
## Docking

To generate in silico structures in KinCo, we started with the DTC dataset (https://drugtargetcommons.fimm.fi/) and identified all the kinase-compound pairs with a binding constant. For each kinase-compound pair, we docked the compound into the homology models of the kinases. The use of homology model allows us to 1) generate in silico structures for the kinase-compound pairs in which the kinase does not yet have a solved structure and 2) sample the conformation state of the kinase. About 11,000 docked poses were generated for each kinase-compound pair in various conformation of the kinases. 

To dock a compound into a kinase, the script dockpipe_qvina.py should be used. 

Information on the coordinates of the poses and the corresponding homology model docked into were saved in the pose.df file which can be loaded with
```python
import pandas as pd
poses = pd.read_pickle(pairbase + 'pose.df')
```
   
   

## Homology Modeling
Homology models used during docking is generated using the Ensembler package https://ensembler.readthedocs.io/en/latest/ by the Chodera lab. This includes the following steps:
1) Gather the sequences of all kinase domains for which we want to generate homology models for (command in ensembler_pipeline.sh)
2) Gather the structures and the sequences of all mammalian kinases which will be used as templates for homology modeling (command in ensembler_pipeline.sh)
3) Align the target kinase and template kinase sequences (command in run_homology_modeling.py)
4) Generate homology models (command in run_homology_modeling.py)
5) Quality check (command in run_homology_modeling.py)

The example resulting files containing the target kinase sequences and the template kinase sequences from step 1 and 2 are included in the target/ and template/ directory respectively.
   
   
   
## Modeling 

### Pretraining on PDBBIND
The model was pretrained on all protein-ligand (not limited to kinases) complexes in PDBBIND2018 dataset.

#### Featurization
Prior to model training, the protein-ligand complexes in PDBBIND2018 in .pdb format were processed (see Processing section below) and featurized using the script feat_PDBBind.py. This script is adapted from https://gitlab.com/cheminfIBB/pafnucy/-/blob/master/prepare.py and uses openbabel to extract atomistic features.

#### Modeling training
To pretrain the models, pretraining_PDBBIND.py was used to train a model using all protein-ligand complexes in PDBBIND2018. 

example training command:
```bash
# point to directory containing the training.hdf, validation.hdf and test.hdf 
input_dir='~/Dataset/PDBBIND2018/cv4'
# point to output directory where the models will be stored under
output_dir='.'

python pretraining_PDBBIND.py --input_dir $input_dir 
			      --output_prefix $output_dir'output' 
   			      --log_dir $output_dir'/logdir/' 
       			      --l2 1e-05 --learning_rate 1e-05 
 			      --rotations 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
			      --num_epochs 20 --num_checkpoints 500 --random_seed 123
```

### Transfer learning using in silico structures
The parameters of the filters were initialized with those learned during pretraining so that the model starts with recognizing biophysically meaningful features derived from crystal structures. 

#### Featurization
To featurize the kinase-compound pairs in PDBBIND2018: feat_xtal.py   
To featurize the kinase-compound pairs selected by vina (M1) or trained models (M2): feat_dock.py
To combine the two datasets: concat_hdfs.py

#### Model training
Models were trained using training_transfer.py

example training command:
```bash
# point to directory containing the training.hdf, validation.hdf and test.hdf               
input_dir='~/Dataset/KinCo/cv4'
# point to output directory where the models will be stored under
output_dir='.'
# whether to freeze the fully connected layer
fine_tune=False
# define the optimal pretrained model to begin training with
meta_path='~/PDBBind2018/cv4/1/output-2020-11-30T15:28:40-232968.meta'

python training_transfer.py --input_dir $input_dir 
			    --output_prefix $output_dir'output' 
			    --log_dir $output_dir'/logdir/' 
			    --l2 1e-05 --learning_rate 1e-05 
			    --num_epochs 20 --num_checkpoints 100 
			    --fine_tune $fine_tune 
			    --random_seed 123 
			    --meta_path $meta_path
```
   
   
  
## Prediction

### Affinity prediction using Autodock Vina
Scoring function from Autodock Vina was used to rescore docked poses produced by Audock Qvina in vina_rescore.py. Affinities produced by Autodock Vina in kcal/mol are converted to Kd in log(uM). This is also the first step to iterative training -- the pose predicted by Vina to have the highest binding affinity was used to train model M1. 
  
### Affinity prediction using trained 3DCNN models
Using trained 3DCNN models to predict binding affinities for a kinase-compound pair. The workflow is consisted of two steps:
1) featurization: feat_poses.py. This script takes multiple kinase-compound pairs and for each pair, featurizes all of its docked poses into a .hdf file for prediction
2) prediction: predict_apair.py. This script takes multiple kinase-compound pairs (listed under input_log/*.txt) and predicts the binding affinities of all the poses using the specified trained models). An example command to run this script is predict.sh

This is the workflow to run the subsequent steps of iterative training in which the representative poses are selected by the models trained in the previous iteration.
   
   
   
## Processing

### process proteins and ligands in PDBBIND2018 for featurization

### distance from center of ATP binding site

### select representative structures for training

#### select pose with highest predicted binding affinity

#### select nonbinding pose as negative examples

