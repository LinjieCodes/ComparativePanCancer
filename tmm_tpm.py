#!/usr/bin/env python

import os
import re
import numpy as np
from rpy2.robjects import r, pandas2ri
r.source('tmm.r')
import pandas as pd

pandas2ri.activate()

#directory of raw count data
count_dir = 'data/count/'

#director to output the TMM data (intermediate file)
tmm_dir = 'data/tmm/'

#director to output the TPM data
tpm_dir = 'data/tpm/'

#read gene exons
exons = {}
with open('gene_annotation/Homo_sapiens.GRCh38.101.gtf') as f:
	for line in f:
		if '\texon\t' in line:
			geneID = re.search('gene_id "(.+?)";', line).group(1)
			if geneID not in exons:
				exons[geneID] = []
			cols = line.split('\t')
			start = int(cols[3])
			end = int(cols[4])
			exons[geneID].append([start, end])
with open('gene_annotation/Mus_musculus.GRCm38.101.gtf') as f:
	for line in f:
		if '\texon\t' in line:
			geneID = re.search('gene_id "(.+?)";', line).group(1)
			if geneID not in exons:
				exons[geneID] = []
			cols = line.split('\t')
			start = int(cols[3])
			end = int(cols[4])
			exons[geneID].append([start, end])
			
#compute gene length
lengths = {}
for geneID in exons:
	sorted_exons = sorted(exons[geneID])
	merged_exon = [sorted_exons[0]]
	for start, end in sorted_exons:
		if start<=merged_exon[-1][1]:
			merged_exon[-1][1] = max(merged_exon[-1][1], end)
		else:
			merged_exon.append([start, end])
	length = 0
	for start, end in merged_exon:
		length += (end-start+1)
	lengths[geneID] = length

#cancer and normal tissues
groups = ['bladder_cancer',
		  'breast_cancer',
		  'colon_cancer',
		  'glioblastoma',
		  'glioma',
		  'kidney_cancer',
		  'leukemia',
		  'liver_cancer',
		  'lung_cancer',
		  'lymphoma',
		  'ovarian_cancer',
		  'pancreatic_cancer',
		  'prostate_cancer',
		  'skin_melanoma',
		  'normal_bladder',
		  'normal_breast',
		  'normal_brain',
		  'normal_colon',
		  'normal_kidney',
		  'normal_blood',
		  'normal_liver',
		  'normal_lung',
		  'normal_ovary',
		  'normal_pancreas',
		  'normal_prostate',
		  'normal_skin']
	
'''
Mouse
'''
#read gene count
counts_list = []
group_samples = {}
for group in groups:
	group_samples[group] = set()
	file = count_dir+'mouse_'+group
	f=open(file)
	counts = {}
	samples = f.readline().strip().split('\t')[1:]
	sampleNum = len(samples)
	for s in samples:
		counts[s] = {}
		group_samples[group].add(s)
	for line in f:
		cols = line.strip().split('\t')
		gene = cols[0]
		exps = [int(e) for e in cols[1:]]
		for i in range(sampleNum):
			s = samples[i]
			counts[s][gene] = exps[i]
	f.close()
	
	#remove sample with poor quality
	counts_poorSampleRemoved = {}
	for s in counts:
		num = 0
		for gene in counts[s]:
			if counts[s][gene] > 10:
				num += 1
		if num > 1000:
			counts_poorSampleRemoved[s] = counts[s]
	
	#convert to dataframe
	counts_df = pd.DataFrame(counts_poorSampleRemoved)	#gene row, sample column
	counts_list.append(counts_df)

#combine into one dataframe
merged_counts_df = pd.concat(counts_list, axis=1) #gene row, sample column
	
#convert to R dataframe
merged_counts_df = pandas2ri.py2rpy(merged_counts_df)

#TMM normalization and write to file
outFile = tmm_dir+'mouse_all_sample'
r.tmm_nor(merged_counts_df, outFile)

#read TMM-normalized counts
tmm_counts = {}
with open(outFile) as f:
	samples = [s.strip('"') for s in f.readline().strip().split('\t')]
	sampleNum = len(samples)
	for s in samples:
		tmm_counts[s] = {}
	for line in f:
		cols = line.strip().split('\t')
		gene = cols[0].strip('"')
		exps = [float(e) for e in cols[1:]]
		for i in range(sampleNum):
			s = samples[i]
			tmm_counts[s][gene] = exps[i]
			
#scale by gene length
lenScaledCounts = {}
for s in tmm_counts:
	lenScaledCounts[s] = {}
	for gene in tmm_counts[s]:
		if gene in lengths:
			lenScaledCounts[s][gene] = tmm_counts[s][gene] / lengths[gene]

#scale by sequencing depth
tpms = {}
for s in lenScaledCounts:
	tpms[s] = {}
	total = sum(lenScaledCounts[s].values())
	for gene in lenScaledCounts[s]:
		tpm = lenScaledCounts[s][gene] * 1000000 / total
		tpms[s][gene] = round(tpm, 2)
			
#output
for group in group_samples:
	outFile = tpm_dir+'mouse_'+group
	this_tpms = {}
	for s in group_samples[group]:
		if s in tpms:
			this_tpms[s] = tpms[s]
	
	this_tpm_df = pd.DataFrame(this_tpms)
	this_tpm_df.to_csv(outFile, sep='\t', index_label='Gene')
	
'''
Human
'''
#read gene count
counts_list = []
group_samples = {}
for group in groups:
	group_samples[group] = set()
	file = count_dir+'human_'+group
	f=open(file)
	counts = {}
	samples = f.readline().strip().split('\t')[1:]
	sampleNum = len(samples)
	for s in samples:
		counts[s] = {}
		group_samples[group].add(s)
	for line in f:
		cols = line.strip().split('\t')
		gene = cols[0]
		exps = [int(e) for e in cols[1:]]
		for i in range(sampleNum):
			s = samples[i]
			counts[s][gene] = exps[i]
	f.close()
	
	#remove sample with poor quality
	counts_poorSampleRemoved = {}
	for s in counts:
		num = 0
		for gene in counts[s]:
			if counts[s][gene] > 10:
				num += 1
		if num > 1000:
			counts_poorSampleRemoved[s] = counts[s]
	
	#convert to dataframe
	counts_df = pd.DataFrame(counts_poorSampleRemoved)	#gene row, sample column
	counts_list.append(counts_df)

#combine into one dataframe
merged_counts_df = pd.concat(counts_list, axis=1) #gene row, sample column
	
#convert to R dataframe
merged_counts_df = pandas2ri.py2rpy(merged_counts_df)

#TMM normalization and write to file
outFile = tmm_dir+'human_all_sample'
r.tmm_nor(merged_counts_df, outFile)

#read TMM-normalized counts
tmm_counts = {}
with open(outFile) as f:
	samples = [s.strip('"') for s in f.readline().strip().split('\t')]
	sampleNum = len(samples)
	for s in samples:
		tmm_counts[s] = {}
	for line in f:
		cols = line.strip().split('\t')
		gene = cols[0].strip('"')
		exps = [float(e) for e in cols[1:]]
		for i in range(sampleNum):
			s = samples[i]
			tmm_counts[s][gene] = exps[i]
			
#scale by gene length
lenScaledCounts = {}
for s in tmm_counts:
	lenScaledCounts[s] = {}
	for gene in tmm_counts[s]:
		if gene in lengths:
			lenScaledCounts[s][gene] = tmm_counts[s][gene] / lengths[gene]

#scale by sequencing depth
tpms = {}
for s in lenScaledCounts:
	tpms[s] = {}
	total = sum(lenScaledCounts[s].values())
	for gene in lenScaledCounts[s]:
		tpm = lenScaledCounts[s][gene] * 1000000 / total
		tpms[s][gene] = round(tpm, 2)
			
#output
for group in group_samples:
	outFile = tpm_dir+'human_'+group
	this_tpms = {}
	for s in group_samples[group]:
		if s in tpms:
			this_tpms[s] = tpms[s]
	
	this_tpm_df = pd.DataFrame(this_tpms)
	this_tpm_df.to_csv(outFile, sep='\t', index_label='Gene')