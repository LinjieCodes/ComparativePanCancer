#!/usr/bin/env python

from scipy import stats

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
		   
for cancer in cancers:
	#read divergent genes
	diverGenes = {}
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			ensg = cols[0]
			ensm = cols[1]
			if abs(float(cols[4]))>1 and float(cols[-1])<0.05:
				diverGenes[ensg] = ensg+'|'+ensm
				diverGenes[ensm] = ensg+'|'+ensm
				
	#read NX in human cancers
	diverGene_NX = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in diverGenes:
				ensg_ensm = diverGenes[cols[0]]
				exps = [float(e) for e in cols[1:]]
				diverGene_NX[ensg_ensm] = exps
	
	#read NX in mouse cancers	
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in diverGenes:
				ensg_ensm = diverGenes[cols[0]]
				exps = [float(e) for e in cols[1:]]
				diverGene_NX[ensg_ensm].extend(exps)
				
	#compute pairwise Spearman correlation
	corrs = []
	finished = set()
	for ensg_ensm1 in diverGene_NX:
		for ensg_ensm2 in diverGene_NX:
			if ensg_ensm1 != ensg_ensm2 and (ensg_ensm1,ensg_ensm2) not in finished:
				exps1 = diverGene_NX[ensg_ensm1]
				exps2 = diverGene_NX[ensg_ensm2]
				corr, p = stats.spearmanr(exps1, exps2)
				corrs.append((ensg_ensm1, ensg_ensm2, corr, p))
				finished.add((ensg_ensm1,ensg_ensm2))
				finished.add((ensg_ensm2,ensg_ensm1))
				
	#write results
	with open('co_expression_network/divergentGene_correlation/'+cancer,'w') as f:
		for ensg_ensm1, ensg_ensm2, corr, p in corrs:
			f.write('\t'.join([ensg_ensm1, ensg_ensm2, str(corr), str(p)])+'\n')