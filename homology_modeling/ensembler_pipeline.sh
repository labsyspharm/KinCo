#!/bin/bash

ensembler init

ensembler gather_targets --gather_from uniprot --query 'kinase AND taxonomy:9606 AND reviewed:yes' --uniprot_domain _regex '^Protein kinase(?!; truncated)(?!; inactive)|PI3K/PI4K|Hexokinase|ADPK|DAGKc|Phosphagen kinase C-terminal|Histidine kinase|Alpha-type protein kinase|RAP'

sbatch -p priority -t 50:00 ensembler gather_templates --gather_from uniprot --query '(family:kinase OR family:hexokinase OR family:glucokinase OR family:PI3/PI4-kinase OR family:"ATP:guanido phosphotransferase" OR family:TAF1) AND taxonomy:40674 AND reviewed:yes' --uniprot_domain _regex '^Protein kinase(?!; truncated)(?!; inactive)|PI3K/PI4K|Hexokinase|ADPK|DAGKc|Phosphagen kinase C-terminal|Histidine kinase|Alpha-type protein kinase|RAP'
