#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:49:40 2020

This script calculates the predicted vina energy. qvina was used for the actual
docking due to its speed, but for energy comparison, we will use the more popular
Autodock vina. For better parallelization, the dataset will be batched. Therefore
the script takes the beginning and the end of the index of the batch as input.
For each kinase-compound pair this scripts calculates the predicted energy vina
predicted for each pose and saves all the poses of this pair in a table

@author: cl321
"""
import numpy as np
import pickle
import pandas as pd
import subprocess
import glob
import os
import sys
import shutil

T = 273+25
R = 8.314
jperkcal = 4184

vina = '/n/groups/alquraishi/Changchang/Softwares/vina/bin/vina'
DeepAlign='/home/cl321/DeepAlign/DeepAlign'
wd = '/n/groups/alquraishi/Changchang/Datasets/symlink_Exp1/'
os.chdir(wd)
wd_tmp = '/n/groups/alquraishi/Changchang/Datasets/vina_energy/tmp/'
biochem_df = pd.read_csv("../resources/DTC_cvsplit1127.csv",
                         low_memory = False, dtype={'ENTREZ_GENE_ID':object,
                                               'INCHIKEY_DESALT':object})
biochem_df.drop_duplicates(subset = ['ENTREZ_GENE_ID', 'INCHIKEY_DESALT'], inplace=True)
biochem_df.reset_index(drop=True, inplace=True)

begin = int(sys.argv[1])
end = int(sys.argv[2])

def unzip_lgd(lgdtmp):
    if os.path.exists(lgdtmp+'.zip'):
        status = subprocess.call(['unzip',lgdtmp+'.zip','-d', lgdtmp])
    if status == 9:
        raise FileNotFoundError
    return lgdtmp+'/'

def kcalpermoltokduM(kcalpermol):
    """ vina outputs energy in kcal/mol, so we need to convert it to log(uM) """
    kdM = np.exp(kcalpermol*jperkcal/(T*R))
    kduM = kdM*1000
    return np.log(kduM)

def synthesize_PDBQT(refpdbqt, coords, outputpdbqt):
    """
    vina only works with PDBQT files, but only the coordinates of the docked
    poses were kept due to size. Here we will regenerate the PDBQT file based 
    on the reference PDBQT ligand file which was used for docking and the 
    coordinates from docking
    input: 
        refpdbqt: the ligand PDBQT file used for docking
        coords: the 3D coordinate of the ligand after docking
        outputpdbqt: desired path of the PDBQT file
    output:
        None, but will directly write the new docked PDBQT file to the output path
        indicated
    """
    # read the reference ligand file to get the atoms table
    with open(refpdbqt, 'r') as f:
        lines = f.readlines()
        f.close()
    
    # plug in the coordinates of the docked pose and save
    acc = 0
    try:
        with open(outputpdbqt, 'w') as PDB:
            for aline in lines:
                if ('ATOM' in aline) or ('UNL' in aline):
                    new_line = list(aline[:])
                    new_coord = [str(i) for i in coords[acc].tolist()]
                    new_line[30:38] = list(str(new_coord[0]).rjust(8))
                    new_line[38:46] = list(str(new_coord[1]).rjust(8))
                    new_line[46:54] = list(str(new_coord[2]).rjust(8))
                    new_str = ''.join(new_line)
                    acc += 1
                else:
                    new_str = aline
                PDB.write(new_str)
    except IndexError:
        raise IndexError

def get_center_coords(geneID, hm_name, inchikey):
    """
    calculate the center of search space for the config file
    same as during docking: the homology model is aligned with a reference kinase
    structure to identify its active site. center of the active site is calculated
    """
    # find the reference target file 1a9u .1
    siteref=('/home/cl321//MoreStructures_v3/1a9u.1.1.1.tgt.pdb')
    recordaln= geneID + '-' + hm_name + '-' + inchikey + '_align'
    outputaln = recordaln+'.pdb'
    
    modeltmp = wd_tmp + geneID + '/PDBaln/' 
    if not os.path.exists(modeltmp):
        os.mkdir(modeltmp)
        
    parsealign= modeltmp + hm_name + '-' + inchikey + '-align_parsed.pdb'
    tgt_location = wd+geneID+'/PDBQT/'+hm_name+'.pdbqt'
    
    if os.path.isfile(parsealign):
        if os.path.getsize(parsealign) < 100:
            os.remove(parsealign)
    
    # align the HM to the reference structure
    if not os.path.isfile(parsealign):
        #runs DeepAlign
        subprocess.call([DeepAlign, (siteref), (tgt_location), "-o", (recordaln)])
        
        #.pdb file actually contains both the pdb of the reference and model, so 
        #here I extract only the structure of the model and put it under its target
        #file
        with open(wd + outputaln, 'r') as f:
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
    
    # from the superposed model, read the center coordinates"""
    with open(parsealign,'r') as g:
        lines=g.readlines()
        for line in lines:
            alist=line.split()
            if alist[1]=="865" and alist[2]=="O" and alist[3]=="MET" and alist[4]=="A":
                try:
                    coord1 = np.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord1=np.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord1=np.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break
        for line in lines:
            alist=line.split()
            if alist[1]=="515" and alist[2]=="NH2" and alist[3]=="ARG" and alist[4]=="A":
                try:
                    coord2 = np.array([float(i) for i in alist[6:9]])
                except ValueError:
                    if len(alist[6])>10:
                        coord2=np.array([float(alist[6][:-8]),float(alist[6][-8:]),float(alist[7])])
                    elif len(alist[7])>10:
                        coord2=np.array([float(alist[6]),float(alist[7][:-8]),float(alist[7][-8:])])
                break        
    center = (coord1+coord2)/2

    return center

def make_config(configpath, geneID, ligand_inchikey, hm_name, pose_ID):
    """
    create the config file for vina. same parameters as during docking was used
    indicates the kinase homology model to run predictions for, center of the 
    search space, size of the search space, and other docking parameters
    """
    # input the homology model to prediction
    tgt_location = wd + geneID + '/PDBQT/' + hm_name + '.pdbqt'
    # calculate the center of search space used during docking -- center of active site
    center = get_center_coords(geneID, hm_name, ligand_inchikey)
    # same search space
    size = 50
    
    with open(configpath, 'w+') as f:
        f.write('receptor='+str(tgt_location)+'\n')
        f.write('  '+'\n')
        f.write('center_x='+str(center[0])+'\n')
        f.write('center_y='+str(center[1])+'\n')
        f.write('center_z='+str(center[2])+'\n')
        f.write('    '+'\n')
        f.write('size_x='+str(size)+'\n')
        f.write('size_y='+str(size)+'\n')
        f.write('size_z='+str(size)+'\n')
        f.write('      '+'\n')
        # same parameters as docking
        f.write('exhaustiveness=8'+'\n')
        f.write('num_modes=100'+'\n')
        f.write('energy_range=3'+'\n')
        f.close()     

def convert_tgt_pdbqt(geneid, hm):
    """regenerate homology model pdbqt if there was an error"""
    if not os.path.exists(wd + geneid + '/PDBQT/error/'):
        os.mkdir(wd + geneid + '/PDBQT/error/')
    print('regenerating homology pdbqt', geneid, hm)
    os.rename(wd + geneid + '/PDBQT/' + hm + '.pdbqt', 
              wd + geneid + '/PDBQT/error/' + hm + '-error.pdbqt')
    
    cmd=["/home/cl321/bin/mglpython",
            "/home/cl321/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py", 
        "-r", wd + geneid + '/PDB/' + hm + '.pdb', 
        '-o', wd + geneid+ '/PDBQT/' + hm + '.pdbqt']
    
    subprocess.call(cmd)
            
def get_kcalpermol(result_string):
    """parse the vina output into kcal/mol
    input: string of the screen output from running vina
    output: energy in kcal/mol
    """
    for aline in result_string.split('\n'):
        if 'Affinity' in aline:
            return float(aline.split()[1])        
    
def calc_affinity(configfile, conformation_path):
    """
    run vina to calculate the affinity based on the parameters in the config file
    input: 
        configfile: path of the configfile
        conformation_path: the pose to predict
    output: 
        energy in kcal/mol and Kd
    """
    # run the vina exe
    run_res = subprocess.run([vina, '--score_only', '--config', configfile,
                              '--ligand', conformation_path], 
                             stdout=subprocess.PIPE)
    # catch the output and parse into kcal and kd
    str_res = run_res.stdout.decode("utf-8") 
    energy_kcal = get_kcalpermol(str_res)
    energy_kd = kcalpermoltokduM(energy_kcal)
    return [energy_kcal, energy_kd]

def process_hm_pose_pair(geneid, inchikey, hm_name, pose_ID):
    """ to predict the energy between a homology model and a pose
    input:
        geneid: gene id of the kinase
        inchikey: inchikey of the compound
        hm_name: homology model of the kinase
        pose_ID: the pose of the ligand 
    output:
        a df with the energy between the homology model and the pose
    """
    # initialize the df via dictionary
    res_dict = {'modelname':[], 'poseID':[], 'energy_kcal':[], 'energy_kd':[]}
    
    # create the config file with parameters for prediction by vina
    configfile = wd_tmp + geneid + '/config/' + hm_name + '-config_qvina.txt'
    if not os.path.exists(configfile):
        make_config(configfile, geneid, inchikey, hm_name, pose_ID)
    
    # point vina to the pose to predict. Because 20 poses exists for each homology
    # model, we will use the same config file for the same homology model and 
    # swap in the path of the poses
    conformation_path = wd_tmp + geneid + '/' + inchikey + '/poses/' + hm_name + '/' + pose_ID + '_qvina_output.pdbqt'
    
    # predict the affinity
    try:
        [energy_kcal, energy_kd] = calc_affinity(configfile, conformation_path)
    except TypeError:
        # error could be due to wrong target pdbqt. regenerate and mark the old pdbqt
        convert_tgt_pdbqt(geneid, hm_name)
        [energy_kcal, energy_kd] = calc_affinity(configfile, conformation_path)
    
    # output results
    res_dict['modelname'].append(hm_name)
    res_dict['poseID'].append(pose_ID)
    res_dict['energy_kcal'].append(energy_kcal)
    res_dict['energy_kd'].append(energy_kd)
    res_df = pd.DataFrame(res_dict)
    return res_df
    
def process_apair(geneid, inchikey):
    """ 
    the main workflow function to calculate the energy of all poses of 
    a kinase-compound pair
    
    input: gene ID of the kinase and inchikey of the compound
    output: the predicted energy by vina for all the poses of the kinase-compound
    pair
    """
    # set up 
    filedir = wd+'/'+geneid+'/'+inchikey 

    if not os.path.exists('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+geneid):
        os.mkdir('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+geneid)
    
    if not os.path.exists(wd_tmp + geneid):
        os.mkdir(wd_tmp + geneid)
    if not os.path.exists(wd_tmp + geneid + '/config'):
        os.mkdir(wd_tmp + geneid + '/config')
        
    if not os.path.exists(wd_tmp + geneid + '/' + inchikey):
        os.mkdir(wd_tmp + geneid + '/' + inchikey)
        os.mkdir(wd_tmp + geneid + '/' + inchikey + '/poses/')
        
    try:
        filedir = unzip_lgd(filedir)
    except FileNotFoundError:
        raise FileNotFoundError

    
    energy_df = pd.DataFrame()
    refpdbqt = filedir + '/PDBQT/' + inchikey + '.pdbqt'
    posedf = pd.read_pickle(filedir+'/pose.df')
    for idx, arow in posedf.iterrows():
        hm = arow['ModelName']
        poseID = arow['PoseID']
        # reproduce the PDBQT file of a pose (one homology model-pose pair)
        outputpdbqt = wd_tmp + geneid + '/' + inchikey + '/poses/' + hm + '/' + poseID + '_qvina_output.pdbqt'
        if not os.path.exists(wd_tmp + geneid + '/' + inchikey + '/poses/' + hm):
            os.mkdir(wd_tmp + geneid + '/' + inchikey + '/poses/' + hm)
        coords = arow['LigandChain'][['x','y','z']].values
        try:
            synthesize_PDBQT(refpdbqt, coords, outputpdbqt)
        except IndexError:
            with open('error_log.txt', 'a') as ef:
                ef.write('{0} {1} {2} {3} pose_error \n'.format(geneid, inchikey, hm, poseID))
            continue
        # calculate the energy of the pose
        pose_energy_df = process_hm_pose_pair(geneid, inchikey, hm, poseID)
        energy_df = pd.concat([energy_df, pose_energy_df])
    
    # save the energy for all the poses
    energy_df.to_csv('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+geneid+'/'+inchikey+'_energy.csv')
    
    if os.path.exists('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+ geneid + '/'+inchikey+'_energy.csv'):
        shutil.rmtree(filedir)
        shutil.rmtree(wd_tmp + geneid + '/' + inchikey)
        
    return energy_df

for idx, arow in biochem_df.iloc[begin:end].iterrows():
    geneid = str(arow['ENTREZ_GENE_ID'])
    lgdinchikey = arow['INCHIKEY_DESALT']
    if not os.path.exists('/n/groups/alquraishi/Changchang/Datasets/vina_energy/'+ geneid + '/'+lgdinchikey+'_energy.csv'):
        try:
            edf = process_apair(geneid, lgdinchikey)
        except FileNotFoundError:
            continue
    
    
    
