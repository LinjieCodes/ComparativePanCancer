#!/usr/bin/env python

from rpy2.robjects import r, pandas2ri
r.source('combat_correct_modeUsed.r')
import os
import pandas as pd
import numpy as np

pandas2ri.activate()

#read human TPM
tpm_dir = 'data/tpm/'
files = [file for file in os.listdir(tpm_dir) if 'human_' in file]
tissue_samples = {}
exps = {}
pre_exps = {}
totalSampleNum = 0
for file in files:
	tissue = file.replace('human_', '')
	with open(tpm_dir+file) as f:
		print(file)
		samples = f.readline().strip().split('\t')[1:]
		tissue_samples[tissue] = samples
		sampleNum = len(samples)
		totalSampleNum += sampleNum
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			if gene not in pre_exps:
				pre_exps[gene] = {}
			for i in range(sampleNum):
				sample = samples[i]
				exp = expList[i]
				pre_exps[gene][sample] = exp
			
#only reserve genes with TPM>0.1 at more than 20% of all samples
for gene in pre_exps:
	if len([e for e in pre_exps[gene].values() if e >0.1]) > 0.2*totalSampleNum:
		for s in pre_exps[gene]:
			if s not in exps:
				exps[s] = {}
			exps[s][gene] = pre_exps[gene][s]
				
#convert to dataframe
#sample columns, gene rows
exps_df = pd.DataFrame(exps)
	
#log2 transform
exps_df = exps_df.applymap(lambda x:np.log2(x+0.01))
	
#all TCGA samples as a batch, all GTEx samples as another
batches = {'studyID':{}}
for sample in exps:
	batches['studyID'][sample] = sample.split('-')[0]
batches_df = pd.DataFrame(batches)	#studyID column, sample row

#group information
groups = {'group':{}}
for group in tissue_samples:
	for sample in tissue_samples[group]:
		groups['group'][sample] = group
groups_df = pd.DataFrame(groups)	#group column, sample row

#call the combat function
cutoff = -10
combat_tpm_dir = 'geneTPM_ComBatCorrected_log2/'
outFile = combat_tpm_dir+'humanAllSample'
result = r.combat_correct(pandas2ri.py2rpy(exps_df), pandas2ri.py2rpy(batches_df), 
						  pandas2ri.py2rpy(groups_df), cutoff, outFile)

#split into files corresponding to tissues
corrected_exps = {}
with open(outFile) as f:
	samples = [s.strip('"') for s in f.readline().strip().split(' ')]
	sampleNum = len(samples)
	for s in samples:
		corrected_exps[s] = {}
	for line in f:
		cols = line.split(' ')
		gene = cols[0].strip('"')
		this_gene_exps = [round(float(e),2) for e in cols[1:]]
		for i in range(sampleNum):
			s = samples[i]
			corrected_exps[s][gene] = this_gene_exps[i]
for tissue in tissue_samples:
	this_tissue_exps = {}
	for s in corrected_exps:
		if s in tissue_samples[tissue]:
			if 'normal' in tissue and 'TCGA-' in s:
				continue
			this_tissue_exps[s] = corrected_exps[s]
			
	#convert to dataframe
	this_tissue_exps_df = pd.DataFrame(this_tissue_exps)
	
	#output
	fileName = 'human_'+tissue
	this_tissue_exps_df.to_csv(combat_tpm_dir+fileName, sep='\t', index_label='Gene')