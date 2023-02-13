#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from sklearn import preprocessing

#read human TPM
print('Reading humam TPM')
human_exps = {}
TPM_dir = 'data/geneTPM_ComBatCorrected_log2/'
files = [file for file in os.listdir(TPM_dir) if 'human_' in file]
for file in files:
	cancer = file.replace('human_', '')
	with open(TPM_dir+file) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleNum = len(samples)
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			for i in range(sampleNum):
				sample = samples[i]+'-'+cancer
				exp = expList[i]
				if sample not in human_exps:
					human_exps[sample] = {}
				human_exps[sample][gene] = exp
			
#read mouse TPM
print('Reading mouse TPM')
mouse_exps = {}
files = [file for file in os.listdir(TPM_dir) if 'mouse' in file]
for file in files:
	cancer = file.replace('mouse_', '')
	with open(TPM_dir+file) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleNum = len(samples)
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			for i in range(sampleNum):
				sample = samples[i]+'-'+cancer
				exp = expList[i]
				if sample not in mouse_exps:
					mouse_exps[sample] = {}
				mouse_exps[sample][gene] = exp
	
#convert to dataframe
#sample row, gene column
print('df convert')
human_exps = pd.DataFrame(human_exps).T
mouse_exps = pd.DataFrame(mouse_exps).T

print('z-score scaling')
z_scaler = preprocessing.StandardScaler()
human_exps = pd.DataFrame(z_scaler.fit_transform(human_exps.values), columns=human_exps.columns, index=human_exps.index)
mouse_exps = pd.DataFrame(z_scaler.fit_transform(mouse_exps.values), columns=mouse_exps.columns, index=mouse_exps.index)
mouse_exps = mouse_exps.dropna(axis=1)	  #clean NaN-containing gene columns

#groups
groups = ['bladder_cancer',
		  'breast_cancer',
		  'glioblastoma',
		  'glioma',
		  'kidney_cancer',
		  'leukemia':',
		  'liver_cancer',
		  'lung_cancer',
		  'lymphoma',
		  'ovarian_cancer',
		  'pancreatic_cancer',
		  'prostate_cancer',
		  'skin_melanoma',
		  'normal_bladder',
		  'normal_blood',
		  'normal_brain',
		  'normal_breast',
		  'normal_kidney',
		  'normal_liver',
		  'normal_lung',
		  'normal_ovary',
		  'normal_pancreas',
		  'normal_prostate',
		  'normal_skin'}

#output the normalized gene expression		  
human_exps_dict = {}
mouse_exps_dict = {}
for group in groups:
	human_exps_dict[group] = {}
	mouse_exps_dict[group] = {}
human_genes = list(human_exps.columns)
mouse_genes = list(mouse_exps.columns)
for sample in human_exps.index:
	this_group = sample.split('-')[-1]
	sampleID = sample.replace('-'+this_group, '')
	exps = list(human_exps.loc[sample])
	human_exps_dict[this_group][sampleID] = {}
	for i in range(len(human_genes)):
		gene = human_genes[i]
		human_exps_dict[this_group][sampleID][gene] = round(exps[i],2)
for sample in mouse_exps.index:
	this_group = sample.split('-')[-1]
	sampleID = sample.replace('-'+this_group, '')
	exps = list(mouse_exps.loc[sample])
	mouse_exps_dict[this_group][sampleID] = {}
	for i in range(len(mouse_genes)):
		gene = mouse_genes[i]
		mouse_exps_dict[this_group][sampleID][gene] = round(exps[i],2)
		
outDir = 'data/scaledGeneExp/'
for group in human_exps_dict:
	file = 'human_'+group
	df = pd.DataFrame(human_exps_dict[group])
	df.to_csv(outDir+file, sep='\t', index_label="Gene")
for group in mouse_exps_dict:
	file = 'mouse_'+group
	df = pd.DataFrame(mouse_exps_dict[group])
	df.to_csv(outDir+file, sep='\t', index_label="Gene")