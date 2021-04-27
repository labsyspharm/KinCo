#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates training/validation/test set partitions for a given list
of kinase-compound pairs based on various kinase AND compound similarity 
thresholds. 

This is achieved by creating the validation/test sets based on clustering:
1) cluster the kinases at specified sequence identity threshold
2) sample a few kinase clusters
3) cluster the compounds at specified tanimoto similarity threshold
4) sample a few compound clusters with compounds interacting with the selected kinases
5) from each of the sampled compound clusters, sample one kinase-compound pair 
    to form the validation/test set
6) remove all the rest of the kinase-compound pairs in the sampled clusters
7) repeat 1)-6) at various kinase and compound thresholds
8) group all pairs in the unsampled clusters into training set
    
"""
import pandas as pd
import pickle
import numpy as np
import scipy
from scipy import spatial
from scipy.cluster.hierarchy import fcluster, linkage
import random
import math
import sys

# parse the input threshold information
kinase_thresholds = []
for i in sys.argv[1].split(','):
    if i!='random':
        kinase_thresholds.append(int(i))
    else:
        kinase_thresholds.append(i)

cmpd_thresholds = []
for i in sys.argv[2].split(','):
    if i!='random':
        cmpd_thresholds.append(float(i))
    else:
        cmpd_thresholds.append('random')
        
num_sel = sys.argv[3]
size_limit = int(sys.argv[4])

output_name = sys.argv[5]
seed = int(output_name.split('-')[-1])

# load other dataset specific information to aid partitioning
fp_dict = pickle.load(open('dataset/morganfp_2_1024_dict_PDBBIND2018+DTC.p', 'rb'))
tgt_dist_dict_by_gene = pickle.load(open('dataset/tgt_dist_dict_by_gene_v4.p', 'rb'))
kinase_bound_per_cmpd = pickle.load(open('dataset/kinase_bound_per_cmpd.p', 'rb'))
cmpd_bound_per_kinase = pickle.load(open('dataset/cmpd_bound_per_kinase.p', 'rb'))

def classify_by_cluster(cluster_indexes, individuals):
    """create the cluster dictionary which maps the index of the cluster (string)
    to its members (list) """
    
    assert len(cluster_indexes)==len(individuals)
     
    # list the indexes of the clusters to partition into cvsets
    original_clusters = sorted(list(set(cluster_indexes)))
    
    # map individuals into their respective clusters
    cluster_w_members = {}
    for acl in original_clusters:
        cluster_w_members[str(acl)] = []
    for acl, amem in zip(cluster_indexes, individuals):
        cluster_w_members[str(acl)].append(amem)
        
    return cluster_w_members

def random_clusters(pairlist, tgt_or_cmpd):
    
    pos_map = {'kinase':0, 'cmpd':1}
    
    num_data = {}
    pseudo_cluseters = {}
    
    for idx, apair in pairlist:
        pseudo_cluseters[str(idx+1)] = apair.split('|')[pos_map[tgt_or_cmpd]]
        num_data[str(idx+1)] = 1
            
    return [pseudo_cluseters, num_data, None]

def cluster_lgd(pairlist, threshold):
    """
    take a list of kinase-compound pairs, cluster them by threshold and 
    return a dictionary: {cluster_index: [compounds in that cluster]}, 
    a list: [number of kinase-compound pairs in each cluster (ordered by index)],
    and a dictionary {compound: index of the cluster it belongs}
    """
    cmpd_list = sorted(list(set([i.split('|')[1] for i in pairlist])))
    print('            clustering %s compounds'%len(cmpd_list))
    
    # create the fingerprint matrix according to the order of cmpd_list
    print('            creating compound clusters')
    fp_mtx = []
    for acmpd in cmpd_list:
        fp_mtx.append(fp_dict[acmpd])
    
    # cluster the cmpd_list so that the distance between the centroids of two clusters is 
    # above threshold
    Ylgd = scipy.spatial.distance.pdist(fp_mtx, metric='jaccard')
    Zlgd = scipy.cluster.hierarchy.linkage(Ylgd, method='single')
        
    cmpds_by_cluster_index = scipy.cluster.hierarchy.fcluster(Zlgd, threshold, criterion='distance')
    cluster_w_cmpds = classify_by_cluster(cmpds_by_cluster_index, cmpd_list)
    print('            produced %s compound clusters'%max(cmpds_by_cluster_index))
    # for each cluster, calculate the number of data points in that cluster (kinases and the associated pairs)
    num_data = {}
    sorted_clusters = range(1, len(cluster_w_cmpds.keys())+1)
    for acl in sorted_clusters:
        data_by_cmpd = [kinase_bound_per_cmpd['num'][i] for i in cluster_w_cmpds[str(acl)]]
        num_data[str(acl)] = sum(data_by_cmpd)
    
    cmpd_to_cmpd_cluster = {}
    for acmpd, cl_idx in zip(cmpd_list, cmpds_by_cluster_index):
        cmpd_to_cmpd_cluster[acmpd] = str(cl_idx)
        
    print('            finish clustering compounds')
    return [cluster_w_cmpds, num_data, cmpd_to_cmpd_cluster]

def tgt_distance(gene1, gene2):
    return tgt_dist_dict_by_gene[frozenset({gene1[0], gene2[0]})]

def cluster_tgt(pairlist, threshold):
    """
    take a list of kinase-compound pairs, cluster them by threshold and 
    return a dictionary: {cluster_index: [kinases in that cluster]} and 
    a list: [number of kinase-compound pairs in each cluster (ordered by index)]
    """
    gene_list = sorted(list(set([i.split('|')[0] for i in pairlist])))
    print('clustering %s proteins'%len(gene_list))
    genes_rs = np.reshape(gene_list, [len(gene_list),1])
        
    # cluster the gene_list so that the distance between the centroids of two clusters is 
    # above threshold
    Ytgt = scipy.spatial.distance.pdist(genes_rs, metric=tgt_distance)
    Ztgt = scipy.cluster.hierarchy.linkage(Ytgt, method='single')
    kinases_by_cluster_index = scipy.cluster.hierarchy.fcluster(Ztgt, threshold, criterion='distance')
    
    cluster_w_kinases = classify_by_cluster(kinases_by_cluster_index, gene_list)
    print('produced %s kinase clusters'%max(kinases_by_cluster_index))
    
    num_data = {}
    sorted_clusters = range(1, len(cluster_w_kinases.keys())+1)
    for acl in sorted_clusters:
        data_by_kinase = [cmpd_bound_per_kinase['num'][i] for i in cluster_w_kinases[str(acl)]]
        num_data[str(acl)] = sum(data_by_kinase)
    
    kinase_to_kinase_cluster = {}
    for acmpd, cl_idx in zip(gene_list, kinases_by_cluster_index):
        kinase_to_kinase_cluster[acmpd] = str(cl_idx)
        
    print('finish clustering proteins')    
    return [cluster_w_kinases, num_data, kinase_to_kinase_cluster]

def random_sample_compound(input_pairs, sampled_kinases, 
                           kinase_to_kinase_cluster, num_sel, kcl_to_remove_acc, ref_pairs):

    num_sel_dict = {'high': int(num_sel.split(',')[0]), 
                    'low': int(num_sel.split(',')[1]), 
                    'xtal':int(num_sel.split(',')[2]),
                    'nonbind': int(num_sel.split(',')[3])}
    
    sampled_pairs_acc = []
    other_cluster_member_record = []
    
    for subset in ['xtal', 'nonbind', 'high', 'low']:
        
        input_pairs_to_sample = sorted(set(ref_pairs[subset])&set(input_pairs))
        
        # identfy the ligands interacting with the kinases of interest. To conserve 
        # data, randomly sample from kinases that have been sampled previously
        interacting_pairs = []
        for akin in sampled_kinases:
            if kinase_to_kinase_cluster[akin] in kcl_to_remove_acc:
                for apair in input_pairs_to_sample:
                    if apair.split('|')[0] == akin:
                        interacting_pairs.append(apair)
                    
        print('            {0} dataset: random sampling from {1} pairs'.format(subset, len(interacting_pairs)))
        
        if len(interacting_pairs)<num_sel_dict[subset]:
            print('            not enough ligands for sampling')
            print('      ')
            if (subset == 'low') or ((num_sel_dict[subset] - len(interacting_pairs)>3)&(len(interacting_pairs)!=0)):
                return [input_pairs, None, None, None]
            else:
                toadd = num_sel_dict[subset] - len(interacting_pairs)
                num_sel_dict['low'] = num_sel_dict['low'] + toadd
                num_sel_dict[subset] = len(interacting_pairs)

        sampled_pairs = random.sample(interacting_pairs, num_sel_dict[subset])
        sampled_pairs_acc += sampled_pairs
        other_cluster_member_record += sampled_pairs
        
        update_input_pairs = sorted(list(set(input_pairs)-set(sampled_pairs_acc)))
            
        kcl_to_remove = []
        for apair in sampled_pairs_acc:
            kcl_to_remove.append(kinase_to_kinase_cluster[apair.split('|')[0]])
        kcl_to_remove = sorted(list(set(kcl_to_remove)))
        
    print('           ',len(sampled_pairs_acc), 'pairs sampled at random')    
    print('            ending with %s pairs'% len(update_input_pairs)) 
    print('      ')
    
    return [update_input_pairs, sampled_pairs_acc, kcl_to_remove, other_cluster_member_record]
    
def sample_compound(input_pairs, cmpd_threshold, 
                    sampled_kinases, kinase_to_kinase_cluster,
                    size_limit, num_sel, ref_pairs):
    
    # cluster the ligands            
    print('            clustering ligands at %s...'%cmpd_threshold)
    print('            starting with %s pairs'% len(input_pairs))      
    
    cmpds_clusters, cmpd_cl_size, cmpd_to_cmpd_cluster = cluster_lgd(input_pairs, cmpd_threshold)
    
    ccl_to_remove = []
    kcl_to_remove = [] # just in case some kinase clusters are not sampled, I want to recycle and use them later
    sampled_pairs_acc = []
    other_cluster_member_record = []
    
    num_sel_dict = {'high': int(num_sel.split(',')[0]), 
                    'low': int(num_sel.split(',')[1]), 
                    'xtal':int(num_sel.split(',')[2]),
                    'nonbind': int(num_sel.split(',')[3])}
    
    # sample public and xtal data separately
    for subset in ['xtal', 'nonbind', 'high', 'low']:
        
        input_pairs_to_sample = sorted(set(ref_pairs[subset])&set(input_pairs))
    
        # identfy the ligands interacting with the kinases of interest
        interacting_ligands = []
        for akin in sampled_kinases:
            for apair in input_pairs_to_sample:
                if apair.split('|')[0] == akin:
                    interacting_ligands.append(apair.split('|')[1])
        interacting_ligands = sorted(list(set(interacting_ligands)))
               
        # Sample the corresponding ligand clusters
            # avoid sampling ligand clusters that are too big to prevent data loss
        interacting_ligand_cluster = sorted(list(set([cmpd_to_cmpd_cluster[algd] for algd in interacting_ligands])))
        print('            {0} dataset: {1} interacting ligands from {2} interacting clusters'.format(subset, len(interacting_ligands), len(interacting_ligand_cluster)))
        
        if len(interacting_ligand_cluster) <= num_sel_dict[subset]:
            print('            not enough ligands for sampling')
            print('      ')
            if (subset == 'low') or ((num_sel_dict[subset] - len(interacting_ligand_cluster)>3)&(len(interacting_ligand_cluster)!=0)):
                return [input_pairs, None, None, None]
            else:
                toadd = num_sel_dict[subset] - len(interacting_ligand_cluster)
                num_sel_dict['low'] = num_sel_dict['low'] + toadd
                num_sel_dict[subset] = len(interacting_ligand_cluster)
                print('            adjusted sampling size: {0}:{1}, {2}:{3}, {4}:{5}, {6}:{7}'.format('xtal', num_sel_dict['xtal'],
                      'high', num_sel_dict['high'],  'low', num_sel_dict['low'], 'nonbind', num_sel_dict['nonbind']))
                if len(interacting_ligand_cluster) == 0:
                    continue
                
        interacting_ligand_cluster_size= {}
        for ccl in interacting_ligand_cluster:
            pair_in_ccl = []
            for acmpd in cmpds_clusters[ccl]:
                pair_in_ccl += kinase_bound_per_cmpd['kinase'][acmpd]
            interacting_ligand_cluster_size[ccl] = len(pair_in_ccl)
            # statistically determine size limit to use
        interacting_ligand_cluster_filtered = []
        stat_cutoff = sorted(interacting_ligand_cluster_size.values())[max(num_sel_dict[subset]-1, round(len(interacting_ligand_cluster)/2))]
        soft_cutoff = min(stat_cutoff, size_limit)
        print('            using size cutoff: %s pairs' %soft_cutoff)
        for ccl, asize in interacting_ligand_cluster_size.items():
            if asize <= soft_cutoff:
                interacting_ligand_cluster_filtered.append(ccl)
            # Do not sample the same cluster in private dataset if already sampled in public
        interacting_ligand_cluster_filtered = sorted(list(set(interacting_ligand_cluster_filtered) - set(ccl_to_remove)))
            # sample the passed ligand clusters
        print('            sampling ligands...')
        print('            %s ligand clusters remainining after filtering'%len(interacting_ligand_cluster_filtered))
        if len(interacting_ligand_cluster_filtered)<num_sel_dict[subset]:
            if (subset == 'low') or (num_sel_dict[subset] - len(interacting_ligand_cluster_filtered)>3):
                print('            not enough ligands for sampling')
                print('      ')
                return [input_pairs, None, None, None]
            else:
                toadd = num_sel_dict[subset] - len(interacting_ligand_cluster_filtered)
                num_sel_dict['low'] = num_sel_dict['low'] + toadd
                num_sel_dict[subset] = len(interacting_ligand_cluster_filtered)
                print('            adjusted sampling size: {0}:{1}, {2}:{3}, {4}:{5}, {6}:{7}'.format('xtal', num_sel_dict['xtal'],
                      'high', num_sel_dict['high'],'low', num_sel_dict['low'], 'nonbind', num_sel_dict['nonbind']))

        sampled_pairs = []
        sampled_cmpd_clusters = random.sample(interacting_ligand_cluster_filtered, num_sel_dict[subset])
        sampled_cmpd_clusters = sorted(list(set(sampled_cmpd_clusters)))
        assert len(set(sampled_cmpd_clusters)) == num_sel_dict[subset]
        ccl_to_remove += sampled_cmpd_clusters
        for accl in sampled_cmpd_clusters:
            # to avoid the performance biases toward clusters with more members, 
            # only sample one pair per ligand cluster
            sampled_cmpds = sorted(list(set(cmpds_clusters[accl])&set(interacting_ligands)))
            random.shuffle(sampled_cmpds)
            random.shuffle(sampled_kinases)
            # flatten ligands into pairs and mark the kinases sampled for later removal
            for acmpd in sampled_cmpds:
                for akin in sampled_kinases:
                    apair = akin+'|'+acmpd
                    if apair in input_pairs_to_sample:
                        sampled_pairs.append(apair)
                        kcl_to_remove.append(kinase_to_kinase_cluster[akin])
                        break
                break    
                
        sampled_pairs = sorted(list(set(sampled_pairs)))
        
        if len(sampled_pairs)<num_sel_dict[subset]:
            print('            not enough ligands sampled')
            print('      ')
            return [input_pairs, None, None, None]
        sampled_pairs_acc += sampled_pairs
        print('           ',len(sampled_pairs), 'pairs sampled at %s'%cmpd_threshold)         
        
        # save other members of the sampled cluster
        for accl in sampled_cmpd_clusters:
            # to avoid the performance biases toward clusters with more members, 
            # only sample one pair per ligand cluster
            sampled_cmpds = sorted(list(set(cmpds_clusters[accl])&set(interacting_ligands)))
            # flatten ligands into pairs and mark the kinases sampled for later removal
            for acmpd in sampled_cmpds:
                for akin in sampled_kinases:
                    apair = akin+'|'+acmpd
                    if apair in input_pairs_to_sample:
                        other_cluster_member_record.append(apair)    
                
    # remove ligand cluster and pairs already sampled from further consideration
    update_input_pairs = []
    for apair in input_pairs:
        [akin, acmpd] = apair.split('|')
        if cmpd_to_cmpd_cluster[acmpd] not in ccl_to_remove:
            update_input_pairs.append(apair)
            
    print('            ending with %s pairs after removing ligands'% len(update_input_pairs))  
    print('   ')
      
    return [update_input_pairs, sampled_pairs_acc, kcl_to_remove, other_cluster_member_record]

def sample_kinase_cmpd_pair(input_pairs, kinase_threshold,
                            sampled_kinase_clusters, 
                            kinase_clusters, kinase_to_kinase_cluster,
                            compound_thresholds, size_limit, num_sel, ref_pairs):
    
    sampled_kinases = []
    for akcl in sampled_kinase_clusters:
        sampled_kinases += kinase_clusters[akcl]
    
    kcl_to_remove_acc = []
    updated_input_pairs = input_pairs[:]
    cmpd_thresh_cv = {}
    
    print('        {0} kinases from {1} sampled kinase clusters'.format(len(set(sampled_kinases)), len(sampled_kinase_clusters)))
    print('        sample clusters: ', ', '.join(sorted(list(sampled_kinase_clusters))[:3]))
    
    for cmpd_threshold in compound_thresholds:
        
        if cmpd_threshold == 'random':
            [res_pairs, sel_pairs, kcl_to_remove, other_cl_mem_record] = random_sample_compound(
                                                        updated_input_pairs, sampled_kinases, 
                                                        kinase_to_kinase_cluster,
                                                        num_sel, kcl_to_remove_acc, ref_pairs)
        else:
            [res_pairs, sel_pairs, kcl_to_remove, other_cl_mem_record] = sample_compound(updated_input_pairs, 
                                            cmpd_threshold, 
                                            sampled_kinases, kinase_to_kinase_cluster,
                                            size_limit, num_sel, ref_pairs)
        updated_input_pairs = res_pairs[:]
        
        cmpd_thresh_cv[str(cmpd_threshold)] = sel_pairs
        if sel_pairs == None:
            return [cmpd_thresh_cv, updated_input_pairs, other_cl_mem_record]
        
        kcl_to_remove_acc = set(kcl_to_remove_acc)|set(kcl_to_remove)    
    
    # update input pairs to remove pairs in the same kinase cluster, skip this
    # step if kinasethreshold is random
    if kinase_threshold == 'random':
        remaining_pairs = updated_input_pairs[:]
    else:
        remaining_pairs = []
        for apair in updated_input_pairs:
            akin = apair.split('|')[0]
            if kinase_to_kinase_cluster[akin] not in kcl_to_remove_acc:
                remaining_pairs.append(apair)
    print('        ending with %s pairs after removing kinases'% len(remaining_pairs)) 
    print('   ')
    
    return [cmpd_thresh_cv, remaining_pairs, other_cl_mem_record]    

def random_kinase_clusters(input_pair):
    kinase_clusters = {}
    kinase_cl_size = {}
    kinase_to_kinase_cluster = {}
    
    for idx, apair in enumerate(input_pair):
        akinase = apair.split('|')[0]
        kinase_clusters[str(idx+1)] = [akinase] 
        kinase_cl_size[str(idx+1)] = 1
        # map the kinase to a non-existing cluster so that it won't be removed
        kinase_to_kinase_cluster[akinase] = str(len(input_pair)*2+idx)
    
    return [kinase_clusters, kinase_cl_size, kinase_to_kinase_cluster]

def sample_for_cv(input_pairs, kinase_threshold, compound_thresholds, 
                            size_limit, num_sel):
    
    try_map = {'PDBBIND':3, 'docked':10}
    random_map = {'PDBBIND':300, 'docked':33}
    
    ref_dict = {'PDBBIND':{'high':[], 
                           'low':pdbbind_pairs, 'nonbind':[], 'xtal':[]}, 
                'docked':{'high':external_high_aff_pairs, 
                          'low':external_low_aff_pairs, 
                          'nonbind': external_nonbinding_pairs, 
                          'xtal':xtal_pairs}}
    
    print('**************************************************')
    print('clustering kinases at %s...'% kinase_threshold)
    if kinase_threshold == 'random':
        kinase_clusters, kinase_cl_size, kinase_to_kinase_cluster = random_kinase_clusters(input_pairs)
    else:
        kinase_clusters, kinase_cl_size, kinase_to_kinase_cluster = cluster_tgt(input_pairs, kinase_threshold)
    kinase_cluster_list = sorted(list(kinase_clusters.keys()))
    
    print('filtering kinase clusters...')
    # PDBBind cluster is generally smaller because kinases are far apart, so random 
    # sampling will sample too many PDBBind. To avoid such, sample PDBBind and 
    # dock poses separately 
    PDBbind_kinase_clusters = []
    docked_kinase_clusters = []
    # classify clusters based on the majority of kinases 
    for akcl in kinase_cluster_list:
        kinase_in_cluster = kinase_clusters[akcl]
        tmp = []
        for akin in kinase_in_cluster:
            if akin in docked_kinases:
                tmp.append(akin)
            else:
                continue
        if len(tmp)/len(kinase_in_cluster) > 0.5:
            docked_kinase_clusters.append(akcl)
        else:
            PDBbind_kinase_clusters.append(akcl)
    print(len(PDBbind_kinase_clusters), 'PDBBind clusters', len(docked_kinase_clusters), 'docked clusters')
    print('         ')
    
    # avoid sampling from kinase clusters too big
    docked_kinase_clusters_filtered = [akcl for akcl in docked_kinase_clusters if kinase_cl_size[akcl]<=size_limit]
    PDBBind_kinase_clusters_filtered = [akcl for akcl in PDBbind_kinase_clusters if kinase_cl_size[akcl]<=size_limit]
    
    kinase_cluster_record = {'docked':{'filtered':docked_kinase_clusters_filtered,
                                       'original':docked_kinase_clusters},
                             'PDBBIND':{'filtered':PDBBind_kinase_clusters_filtered,
                                        'original':PDBbind_kinase_clusters}}
    res_dict = {}
    record_dict = {}
    
    remaining_pairs = input_pairs[:]
    # sample 10% of kinase clusters
    for atype in ['docked', 'PDBBIND']:
        print('    working on %s...' %atype)
        res_dict[atype] = {'validation':None, 'test':None}
        record_dict[atype] = {'validation':None, 'test':None}
        
        # checking if we have enough data for cv set before proceeding
        original_length = len(kinase_cluster_record[atype]['original'])
        if kinase_threshold == 'random':
            num_kcl_to_sample = random_map[atype]
        else:
            num_kcl_to_sample = round(0.1*original_length)
        if len(kinase_cluster_record[atype]['filtered']) < 2*num_kcl_to_sample:
            print('not enough kinase clusters after filtering', kinase_threshold)
            return [res_dict, input_pairs, None] 

        print('        making %s validation set'%atype)
        print('     ')
        # when ligand sampling fails, try another group of kinases
        for atry in range(try_map[atype]):
            print('        running try', atry)
            kinase_cluster_record[atype]['validation'] = random.sample(kinase_cluster_record[atype]['filtered'], num_kcl_to_sample)
            [val_cmpd_thresh_cv, res_pairs, val_other_pair_record] = sample_kinase_cmpd_pair(remaining_pairs, 
                                                kinase_threshold,
                                                kinase_cluster_record[atype]['validation'],
                                                kinase_clusters, kinase_to_kinase_cluster,
                                                compound_thresholds, size_limit, num_sel, ref_dict[atype])
            if None not in val_cmpd_thresh_cv.values():
                break
        
        remaining_pairs = res_pairs[:]    
        res_dict[atype]['validation']=val_cmpd_thresh_cv
        record_dict[atype]['validation']=val_other_pair_record
        if None in val_cmpd_thresh_cv.values():
            print('        not enough pairs for validation set')
            break
        
        print('        making %s test set '%atype)
        print('     ')
        # when ligand sampling fails, try another group of kinases
        for atry in range(try_map[atype]):
            print('        running try', atry)
            kinase_cluster_list_filtered_remaining = sorted(list(set(kinase_cluster_record[atype]['filtered'])-set(kinase_cluster_record[atype]['validation'])))
            kinase_cluster_record[atype]['test'] = random.sample(kinase_cluster_list_filtered_remaining, num_kcl_to_sample)
            [test_cmpd_thresh_cv, res_pairs, test_other_pair_record] = sample_kinase_cmpd_pair(remaining_pairs, 
                                                    kinase_threshold,
                                                    kinase_cluster_record[atype]['test'], 
                                                    kinase_clusters, kinase_to_kinase_cluster,
                                                    compound_thresholds, size_limit, num_sel, ref_dict[atype])
            if None not in test_cmpd_thresh_cv.values():
                break
            
        remaining_pairs = res_pairs[:]    
        res_dict[atype]['test']=test_cmpd_thresh_cv
        record_dict[atype]['test']=test_other_pair_record
        if None in test_cmpd_thresh_cv.values():
            print('        not enough pairs for test set')
            break
        
    return [res_dict, remaining_pairs, record_dict]

def create_cv(input_pairs, kinase_thresholds, compound_thresholds, size_limit, num_sel):
    
    remaining_pairs = input_pairs[:]
    
    kinase_thresh_cv = {}
    cvrecord_kinase_thresh = {}
    
    for akthresh in kinase_thresholds:
        kinase_thresh_cv[str(akthresh)] = {}
        cvrecord_kinase_thresh[str(akthresh)] = {}
        
        [val_test_cvdict, remaining_pairs, val_test_record_dict] = sample_for_cv(remaining_pairs, akthresh, compound_thresholds, 
                                                            size_limit, num_sel)

        for acv in ['validation', 'test']:
            kinase_thresh_cv[str(akthresh)][acv] = {}
            cvrecord_kinase_thresh[str(akthresh)][acv] = {}
            
            for atype in ['docked', 'PDBBIND']:   
                try:
                    kinase_thresh_cv[str(akthresh)][acv][atype] = val_test_cvdict[atype][acv]
                    cvrecord_kinase_thresh[str(akthresh)][acv][atype] = val_test_record_dict[atype][acv]
                except KeyError:
                    kinase_thresh_cv[str(akthresh)][acv][atype] = None
                    cvrecord_kinase_thresh[str(akthresh)][acv][atype] = None
            if None in val_test_cvdict['docked'][acv].values():
                return [kinase_thresh_cv, remaining_pairs, cvrecord_kinase_thresh]
    
    kinase_thresh_cv['training'] = remaining_pairs
    print ('training set size: ', len(kinase_thresh_cv['training']))    
    return [kinase_thresh_cv, cvrecord_kinase_thresh]

def read_pair_strings(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    pairs = [i[:-1] for i in lines]
    return pairs

if __name__ == '__main__':
    
    print('using random seed', seed)
    random.seed(seed)

    # use new function to read files instead of numpy because it cuts off the '-N' in INCHIKEYS
    docked_kinases = read_pair_strings('dataset/docked_kinases.txt')
    pair_list = read_pair_strings('dataset/PDBBIND+DTC_pairs.txt')
        
    external_high_aff_pairs = read_pair_strings('dataset/highaff_pairs.txt')
    external_low_aff_pairs = read_pair_strings('dataset/lowaff_pairs.txt')
    external_nonbinding_pairs = read_pair_strings('dataset/nonbinding_pairs.txt')
    xtal_pairs = read_pair_strings('dataset/xtal-kinase_pairs.txt')
    pdbbind_pairs = read_pair_strings('dataset/xtal-nonkinase_pairs.txt')
            
    cvsplit = create_cv(pair_list, kinase_thresholds, cmpd_thresholds, 
                                size_limit = size_limit, num_sel = num_sel)
    
    if len(cvsplit) == 2:
        [cvsplit_cvset, cvsplit_member_record] = cvsplit
        with open('cvsplit/cvsplit_%s.p'%output_name, 'wb') as handle:
            pickle.dump(cvsplit_cvset, handle)
        with open('cvsplit/cvrecord_%s.p'%output_name, 'wb') as handle:
            pickle.dump(cvsplit_member_record, handle)
    elif len(cvsplit) == 3:
        with open('cvsplit/failed_cvsplit_%s.p'%output_name, 'wb') as handle:
            pickle.dump(cvsplit, handle)
        
