#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:25:30 2017

@author: cl321
"""
""" 

This is the pipeline to dock a ligand into a kinase. It should be 
called by the run script run_dock_qvina.py which runs the docking of a pdb structure 
by dockingthe ligand into ALL the homology models of the kinase. As a result, 
this script takes in as argument the name of the ligand, the name of the kinase, 
Uniprot name of the kinase domain, and the number of runs
to dock. 

 """

import glob
import subprocess
import os
import argparse
import shutil
import numpy
import gzip
import pickle
import time
import sys
import pandas as pd
from biopandas.pdb import PandasPdb
from biopandas.mol2 import PandasMol2
from collections import Counter

print (sys.version_info)

parser = argparse.ArgumentParser()
parser.add_argument('ligand_inchikey', help='ligand')
parser.add_argument('geneID', help='kinase gene nmae')
parser.add_argument('subdomain', help='kinase Uniprot domain')
parser.add_argument('runs', help='number of runs')
parser.add_argument('ncpus', help='number of CPUs')
args = parser.parse_args()

""" STEP 0. create directory structures """ 
# set up directories
parentfldr = '/n/groups/alquraishi/Changchang/Ensembler/symlink_Exp1/'
parenttmpfldr = '/n/scratch3/users/c/cl321/Ensembler/symlink_Exp1/'
datasetfldr = '/n/groups/alquraishi/Changchang/Datasets/symlink_Exp1/'
os.chdir(parenttmpfldr)
tgtfldr = parentfldr + '/' + args.geneID + '/'
tgttmpfldr = parenttmpfldr + '/' + args.geneID + '/'
lgdfldr = tgttmpfldr + '/' + args.ligand_inchikey

if not os.path.exists(lgdfldr+'/poses/'):
    os.mkdir(lgdfldr+'/poses/')
    
# softwares to use
obabel = '/n/groups/alquraishi/Changchang/openbabel-2a/bin/obabel' 
fconv = '/n/groups/alquraishi/Apps/fconv_124_bin_linux32/fconv'
DeepAlign='/home/cl321/DeepAlign/DeepAlign'
docking='/n/groups/alquraishi/Changchang/qvina/qvina02'

# load the homology models to run and the number of dockings to be run for each
with open((lgdfldr+'/run_number='+str(int(args.runs))), 'r+') as f1:
    lines = f1.readlines()
    undermodellist = [aline.strip() for aline in lines]
with open((lgdfldr+'/run_number='+str(int(args.runs) + 1)), 'r+') as f2:
    lines = f2.readlines()
    overmodellist = [aline.strip() for aline in lines]
    
""" STEP 1. convert ligand smiles to mol2 and PDBQT files """
if not os.path.isfile((lgdfldr +'/mol2/'+args.ligand_inchikey +'.mol2')):
    with open((lgdfldr+'/smiles'), 'r+') as smi:
        ligand_smiles=smi.readline().strip('\n')
    subprocess.call([obabel, ('-:'+ligand_smiles), '-omol2', '-O', 
                        (lgdfldr +'/mol2/' + args.ligand_inchikey + '.mol2'), '--gen3d'])

lgdpdbqt = lgdfldr +'/PDBQT/'+args.ligand_inchikey +'.pdbqt'
if not os.path.isfile(lgdpdbqt):   
    cmd=["/home/cl321/bin/mglpython",
    "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py", 
    "-l", (lgdfldr +'/mol2/' + args.ligand_inchikey +'.mol2'), '-o', lgdpdbqt]
    subprocess.call(cmd)
    
def Convert_LigandChains(ligand_mol2):
    ligand = PandasMol2().read_mol2(ligand_mol2)
    ligandchain=pd.concat([ligand.df['atom_type'],ligand.df['x'],ligand.df['y'],
                           ligand.df['z']],axis=1)	
    return ligandchain

def SplitPDBQTtoMol2(afile, modeldir):
    subprocess.call([obabel, 
                 "-ipdbqt", afile, "-omol2", "-O" + modeldir + poseindex 
                 + "single.mol2", "-m"])
    conflist = sorted(glob.glob(lgdfldr+'/poses/'+model_name+'/'+poseindex+'single*.mol2'))
    return conflist

def Mol2toMol2(confindex):
    tmpmol2 = lgdfldr + '/'+ poseindex + '-' + confindex + '_m0.mol2'
    subprocess.call([fconv, '-T','--m=0', aconf, '--t='+tmpmol2])
    return tmpmol2

def Docking(model, to_run):
    """main script for docking. steps 2-5
    step 2: locate homology models and kept record
    step 3: prepare homology models and convert into PDBQT for docking
    step 4: identify active site for docking via structural alignment
            generate config file with parameters for docking
    step 5: perform docking based on the config file
    input: homology model to use, number of times to perform docking
    output: (for record keeping) homology model docked, 
            number of times to perform docking
    """
    modeltmp = lgdfldr + '/poses/' + model + '/'
    if not os.path.exists(modeltmp):
        os.mkdir(modeltmp)
        
    """ STEP 2. pick out the model file to dock into """
    # first check if the models have been picked. If so, skip this step
    modelname=(tgtfldr +'/PDB/'+model+'.pdb')
    if not os.path.isfile(modelname):
        modelpdb=('/n/groups/alquraishi/Changchang/Ensembler/database_it2/models/'+args.subdomain+
                    '/'+model+'/model.pdb')
        modelpdbgz=('/n/groups/alquraishi/Changchang/Ensembler/database_it2/models/'+args.subdomain+
                    '/'+model+'/model.pdb.gz')
        try:
            assert os.path.isfile(modelpdb)
        except AssertionError:
            # Write uncompressed model.pdb files from model.pdb.gz if necessary
            with gzip.open(modelpdbgz) as model_pdbfile_compressed:
                with open(modelpdb, 'w') as model_pdbfile:
                    model_pdbfile.write(model_pdbfile_compressed.read())    
        shutil.copy(modelpdb, modelname)
        
    print ('finished picking out model file')
    
    """ STEP 3. convert target PDB file to PDBQT """
    # The PDBQT file will first be generated under /$target/PDB/PDBQT/ directory
    # first check if the pdbqt have been generated. If so, skip this step
    modelpdbqt = (tgtfldr +'/PDBQT/'+ model+'.pdbqt')
    if not (os.path.isfile(modelpdbqt) and os.path.getsize(modelpdbqt)!=0):
        cmd=["/home/cl321/bin/mglpython",
             "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py", 
        "-r",modelname, '-o', modelpdbqt]
        subprocess.call(cmd)

    print ('finished converting to pdb')
    
    """ STEP 4. generate config script """    

    """ STEP 4.0 superpose the model to the reference PDB file (to define binding site) """
    # find the reference target file 1a9u .1
    siteref=('./1a9u.1.1.1.tgt.pdb')
    # align homology model to the reference file
    recordaln= args.geneID + '-' + args.ligand_inchikey + '-' + model+'_align'
    outputaln = recordaln+'.pdb'
    parsealign= modeltmp +'/align_parsed.pdb'
    tgt_location = modelpdbqt
    if (not os.path.isfile(parsealign)) or (os.path.getsize(parsealign)==0):
        subprocess.call([DeepAlign, (siteref), (tgt_location), "-o", (recordaln)])
        # .pdb file actually contains both the pdb of the reference and model, so 
        # here I extract only the structure of the model and put it under its target
        # file
        with open((parenttmpfldr+outputaln)) as f:
            lines=f.readlines()
            with open(parsealign, 'a+') as g:
                for line in lines:
                    if 'TER' in line:
                        break
                    g.write(line)
        #delete the temporary files
        outputalnlist=[(recordaln+'.fasta'), (recordaln+'.local'),(recordaln+'.score'),(recordaln+'.pdb')]
        for anoutput in outputalnlist:    
            os.remove(anoutput)
    
    """ STEP 4.1 from the superposed model, read the center coordinates"""
    # binding site is defined by the R515 and M865 in 1A9U
    with open(parsealign,'r') as g:
        lines=g.readlines()
        for line in lines:
            alist=line.split()
            if alist[1]=="865" and alist[2]=="O" and alist[3]=="MET" and alist[4]=="A":
                try:
                    coord1 = numpy.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord1=numpy.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord1=numpy.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break
        for line in lines:
            alist=line.split()
            if alist[1]=="515" and alist[2]=="NH2" and alist[3]=="ARG" and alist[4]=="A":
                try:
                    coord2 = numpy.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord2=numpy.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord2=numpy.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break        
    center = (coord1+coord2)/2
    
    os.remove(parsealign)
        
    """ STEP 4.2 write the configuration files for docking """            
    lgd_location=lgdfldr + '/PDBQT/'+ args.ligand_inchikey + '.pdbqt'
    
    # TODO: search space was intentially kept large to produce negative poses and for 
    # allosteric inhibitors. In the future can narrow down to focus on ATP-competitive
    # inhibitors binding to the active site (eg. size = 10)
    size=50 
    
    configfile=(modeltmp + '/config_qvina'+'.txt')
    if not os.path.exists(configfile):
        f=open(configfile, 'w')
        f.write('receptor='+str(tgt_location)+'\n')
        f.write('ligand='+str(lgd_location)+'\n')
        f.write('  '+'\n')
        f.write('center_x='+str(center[0])+'\n')
        f.write('center_y='+str(center[1])+'\n')
        f.write('center_z='+str(center[2])+'\n')
        f.write('    '+'\n')
        f.write('size_x='+str(size)+'\n')
        f.write('size_y='+str(size)+'\n')
        f.write('size_z='+str(size)+'\n')
        f.write('      '+'\n')
        f.write('exhaustiveness=8'+'\n')
        f.write('num_modes=100'+'\n')
        f.write('energy_range=3'+'\n')
        f.close()
    
    print ('finished writing config files')
        
    """ STEP 5. docking """    
    for i in range(to_run):
        cmd = [docking, '--cpu', args.ncpus, '--config', configfile, '--out', (modeltmp 
             + str(i)+'_qvina_output.pdbqt')]
        subprocess.call(cmd)
        
    return [model, to_run]

# actual step to do the dockng
for amodel in undermodellist: 
    runs = int(args.runs)
    checklist = Docking(amodel, runs)
    print ('finished ', checklist[1], 'runs ', checklist[0])
if overmodellist != []:
    for amodel in overmodellist: 
        runs = int(args.runs) + 1
        checklist = Docking(amodel, runs)
        print ('finished ', checklist[1], 'runs ', checklist[0])

# process docked poses into ligand chain
filelist = sorted(glob.glob(lgdfldr + '/poses' + '/*/*_qvina_output.pdbqt'))
ran = len(filelist)

    # raise warning if number of poses doesn't match 
if ran<560:
    print ('Warning: not enough poses', len(filelist))
    with open(lgdfldr+'/fewer', 'w+') as w:
        w.write(str(len(filelist)))
elif ran>560:
    print ('Warning: too many poses', len(filelist))
    with open(lgdfldr+'/more', 'w+') as w:
        w.write(str(len(filelist)))
else:
    if os.path.exists(lgdfldr+'/fewer'):
        os.remove(lgdfldr+'/fewer')
    elif os.path.exists(lgdfldr+'/more'):
        os.remove(lgdfldr+'/more')

# write all ligand chains into a dataframe
# Atom typing scheme from fconv will be used, so need to convert the files to 
# fconv mol2 format
# !!have to first convert PDBQT to obabel mol2 then fconv mol2, otherwise
# the order of the atoms will not be conserved
d = {'ModelName':[], 'PoseID':[], 'LigandChain':[]}
d_df = pd.DataFrame(d)

for afile in filelist:
    model_name = afile.split('/')[-2]
    poseindex = afile.split('/')[-1].split('_')[0]
    modeldir = lgdfldr + '/poses/' + model_name + '/'
    conflist = SplitPDBQTtoMol2(afile, modeldir)
    for aconf in conflist:
        confindex = aconf.split('.')[0].split('single')[1]
        tmpmol2 = Mol2toMol2(confindex)
        ligand_chain = Convert_LigandChains(tmpmol2)
        d_df=pd.concat([d_df, pd.DataFrame({'ModelName':[model_name], 
                                            'PoseID':[poseindex+'-'+confindex], 
                                            'LigandChain':[ligand_chain]})], 
                        ignore_index = True)
        os.remove(aconf)
        os.remove(tmpmol2)
    os.remove(afile)
# clean up and delete intermediate files
shutil.rmtree(lgdfldr + '/poses/')

# save pose.df and .zip
d_df.to_pickle(lgdfldr + '/pose.df')
shutil.make_archive(lgdfldr, 'zip', lgdfldr)

# quality check and move to dataset directory for long term storage
if (os.path.exists(lgdfldr+'.zip')) and (os.path.getsize(lgdfldr+'.zip') > 2000):
    datasettgt = datasetfldr + args.geneID + '/'
    if not os.path.exists(datasettgt):
        os.mkdir(datasettgt)
    print('copying to dataset')
    subprocess.call(['rsync', '-az', lgdfldr + '.zip', datasettgt])
    shutil.rmtree(lgdfldr)
    print('removed ', lgdfldr)
    if os.path.exists(lgdfldr):
        print ('still here ', lgdfldr)
        shutil.rmtree(lgdfldr)
    if os.path.exists(datasettgt+args.ligand_inchikey+'.zip'):
        print ('removing ', lgdfldr + '.zip')
        os.remove(lgdfldr + '.zip')
else:
    print ('zip failed')

