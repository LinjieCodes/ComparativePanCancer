#!/usr/bin/env python

from rpy2.robjects import r, pandas2ri
r.source('combat_correct.r')
import os
import pandas as pd
import numpy as np

pandas2ri.activate()

#read all samples' studyIDs
studies = {}
with open('sample_list_outlier_unfiltered') as f:
	f.readline()
	for line in f:
		study, sample = line.split('\t')[:2]
		studies[sample] = study

#correct mouse gene TPM
tpm_dir = 'data/tpm/'
combat_tpm_dir = 'geneTPM_ComBatCorrected_log2/'
files = [file for file in os.listdir(tpm_dir) if 'mouse_' in file]
for file in files:
	#read gene TPM
	exps = {}
	with open(tpm_dir+file) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleNum = len(samples)
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			for i in range(sampleNum):
				sample = samples[i]
				exp = expList[i]
				if sample not in exps:
					exps[sample] = {}
				exps[sample][gene] = exp
		
	#convert to dataframe
	#sample columns, gene rows
	exps_df = pd.DataFrame(exps)
	
	#log10-transform
	exps_df = exps_df.applymap(lambda x:np.log10(x+0.01))
	
	#study information
	this_studies = {'studyID':{}}
	for sample in exps:
		this_studies['studyID'][sample] = studies[sample]
	this_studies_df = pd.DataFrame(this_studies)	#studyID column, sample row
	
	cutoff = -1.8
	outFile = combat_tpm_dir+file
	result = r.combat_correct(pandas2ri.py2rpy(exps_df), pandas2ri.py2rpy(this_studies_df), cutoff, outFile)