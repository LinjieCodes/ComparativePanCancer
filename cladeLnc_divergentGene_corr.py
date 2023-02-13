#!/usr/bin/env python

from scipy import stats
import re
import numpy as np
import os

#13 cancer types
cancers = ['bladder_cancer',
		   'breast_cancer',
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
		   'skin_melanoma']
		   
#read TPM expression in cancers
#retain lncRNAs with TPM>0.1 in 
#>=50% of samples of each cancer type
humanLncs_expressed = {}
mouseLncs_expressed = {}
for cancer in cancers:
	#extract clade-specific lncRNAs in human cancer
	humanLncs = set()
	with open('lncrna_analysis/cladeLnc_inCancers/human_'+cancer) as f:
		f.readline()
		for line in f:
			humanLncs.add(line.split('\t')[0])
			
	#extract clade-specific lncRNAs in mouse cancer
	mouseLncs = set()
	with open('lncrna_analysis/cladeLnc_inCancers/mouse_'+cancer) as f:
		f.readline()
		for line in f:
			mouseLncs.add(line.split('\t')[0])
					
	#extract divergent genes
	diverGenes_human = set()
	diverGenes_mouse = set()
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			ensg = cols[0]
			ensm = cols[1]
			deltaNX = float(cols[4])
			fdr = float(cols[-1])
			if abs(deltaNX)>1 and fdr<0.05:
				diverGenes_human.add(ensg)
				diverGenes_mouse.add(ensm)
				
	#extract normalized expression in cancers
	humanExps = {}
	mouseExps = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			exps = [float(e) for e in cols[1:]]
			humanExps[cols[0]] = exps
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			exps = [float(e) for e in cols[1:]]
			mouseExps[cols[0]] = exps
			
	#compute Spearman correlation of lncRNA-divergentGene
	with open('co_expression_network/cladeLnc_divergentGene_correlation/human_'+cancer, 'w') as f:
		#human, lncRNA-divergentGene
		corrs = []
		for lnc in humanLncs:
			if lnc not in humanExps:
				continue
			lncExps = humanExps[lnc]
			if len(set(lncExps))==1:
				continue
			for diverGene in diverGenes_human:
				diverGeneExps = humanExps[diverGene]
				if len(set(diverGeneExps))==1:
					continue
				corr, p = stats.spearmanr(lncExps, diverGeneExps)
				corrs.append([lnc, diverGene, corr, p])
		
		'''		
		#FDR correction
		pValues = [a[3] for a in corrs]
		pValues = np.asfarray(pValues)
		by_descend = pValues.argsort()[::-1]
		by_orig = by_descend.argsort()
		steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
		fdrs = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
		for i in range(len(fdrs[by_orig])):
			corrs[i].append(fdrs[by_orig][i])'''
			
		for lnc, diverGene, corr, p in corrs:
			f.write('\t'.join([lnc, diverGene, str(corr), str(p)])+'\n')
				
	with open('co_expression_network/cladeLnc_divergentGene_correlation/mouse_'+cancer, 'w') as f:
		#mouse, lncRNA-divergentGene
		corrs = []
		for lnc in mouseLncs:
			if lnc not in mouseExps:
				continue
			lncExps = mouseExps[lnc]
			if len(set(lncExps))==1:
				continue
			for diverGene in diverGenes_mouse:
				diverGeneExps = mouseExps[diverGene]
				if len(set(diverGeneExps))==1:
					continue
				corr, p = stats.spearmanr(lncExps, diverGeneExps)
				corrs.append([lnc, diverGene, corr, p])
		
		'''
		#FDR correction
		pValues = [a[3] for a in corrs]
		pValues = np.asfarray(pValues)
		by_descend = pValues.argsort()[::-1]
		by_orig = by_descend.argsort()
		steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
		fdrs = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
		for i in range(len(fdrs[by_orig])):
			corrs[i].append(fdrs[by_orig][i])'''
			
		for lnc, diverGene, corr, p in corrs:
			f.write('\t'.join([lnc, diverGene, str(corr), str(p)])+'\n')