#!/usr/bin/env python

from scipy import stats
import numpy as np
import re

#read one-to-one orthologs
orthologs = {}
with open('gene_annotation/one-to-one-orthologs') as f:
	for line in f:
		ensg, ensm = line.strip().split('\t')
		orthologs[ensg] = ensm
		
#read gene symbols
geneSymbols = {}
with open('gene_annotation/Homo_sapiens.GRCh38.101.gtf') as f:
	for line in f:
		if '\tgene\t' in line:
			geneID = re.search('gene_id "(.+?)";', line)
			symbol = re.search('gene_name "(.+?)";', line)
			if geneID and symbol:
				geneID = geneID.group(1)
				symbol = symbol.group(1)
				geneSymbols[geneID] = symbol
with open('gene_annotation/Mus_musculus.GRCm38.101.gtf') as f:
	for line in f:
		if '\tgene\t' in line:
			geneID = re.search('gene_id "(.+?)";', line)
			symbol = re.search('gene_name "(.+?)";', line)
			if geneID and symbol:
				geneID = geneID.group(1)
				symbol = symbol.group(1)
				geneSymbols[geneID] = symbol

#13 cancer types
cancers = {'bladder_cancer':'normal_bladder',
		   'breast_cancer':'normal_breast',
		   'glioblastoma':'normal_brain',
		   'glioma':'normal_brain',
		   'kidney_cancer':'normal_kidney',
		   'leukemia':'normal_blood',
		   'liver_cancer':'normal_liver',
		   'lung_cancer':'normal_lung',
		   'lymphoma':'normal_blood',
		   'ovarian_cancer':'normal_ovary',
		   'pancreatic_cancer':'normal_pancreas',
		   'prostate_cancer':'normal_prostate',
		   'skin_melanoma':'normal_skin'}

#director of the TPM data
tpm_dir = 'data/geneTPM_ComBatCorrected_log2/'

#director of the NX data
nx_dir = 'data/scaledGeneExp/'
		   
for cancer in cancers:
	normal = cancers[cancer]
	
	#extract DEGs in human cancer
	human_degs = set()
	with open('DEG/human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			if abs(float(cols[1]))>1 and float(cols[-1])<0.05:
				human_degs.add(gene)
	
	#compute human median TPM
	human_tpms = {}
	with open(tpm_dir+'human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms), 2)
			human_tpms[gene] = [medianTpm]
	with open(tpm_dir+'human_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms), 2)
			human_tpms[gene].append(medianTpm)
			
	#compute mouse median TPM
	mouse_tpms = {}
	with open(tpm_dir+'mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms), 2)
			mouse_tpms[gene] = [medianTpm]
	with open(tpm_dir+'mouse_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms), 2)
			mouse_tpms[gene].append(medianTpm)
	
	#read human normalized expression
	human_exps = {}
	with open(nx_dir+'human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			human_exps[gene] = exps
			
	#read mouse normalized expression
	mouse_exps = {}
	with open(nx_dir+'mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			mouse_exps[gene] = exps
			
	#two-sided Mann–Whitney test
	betweenSpeciesDegs = []
	for ensg in human_degs:
		if ensg in orthologs:
			ensm = orthologs[ensg]
			if ensm in mouse_exps:
				medianExp1 = np.median(human_exps[ensg])
				medianExp2 = np.median(mouse_exps[ensm])
				diffExp = round(medianExp1-medianExp2, 2)
				p = stats.mannwhitneyu(human_exps[ensg], mouse_exps[ensm])[1]
				betweenSpeciesDegs.append([p, diffExp, ensg, ensm])
		
	#Benjamini-Hochberg p-value correction
	pValues = [a[0] for a in betweenSpeciesDegs]
	pValues = np.asfarray(pValues)
	by_descend = pValues.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
	fdrs = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
	for i in range(len(fdrs[by_orig])):
		betweenSpeciesDegs[i].append(fdrs[by_orig][i])
		
	#output interspecies DEG
	with open('divergence_result/'+cancer, 'w') as f:
		f.write('GeneID\tGeneID\tGeneSymbol\tGeneSymbol\tdeltaNX\tCancerlog2TPM_human\tNormallog2TPM_human\tCancerlog2TPM_mouse\tNormallog2TPM_mouse\tPValue\tFDR\n')
		for p, diffExp, ensg, ensm, fdr in betweenSpeciesDegs:
			cancermedianTPM_human, normalmedianTPM_human = human_tpms[ensg]
			cancermedianTPM_mouse, normalmedianTPM_mouse = mouse_tpms[ensm]
			geneSymbol_human = geneSymbols[ensg]
			geneSymbol_mouse = geneSymbols[ensm]
			f.write('\t'.join([ensg, ensm, geneSymbol_human, geneSymbol_mouse, 
					str(diffExp), str(cancermedianTPM_human), 
					str(normalmedianTPM_human), str(cancermedianTPM_mouse), 
					str(normalmedianTPM_mouse), str(p), str(fdr)])+'\n')