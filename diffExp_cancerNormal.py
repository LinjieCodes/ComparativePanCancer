#!/usr/bin/env python

from scipy import stats
import numpy as np

#director of the TPM data
tpm_dir = 'data/geneTPM_ComBatCorrected_log2/'

#director of the NX data
nx_dir = 'data/scaledGeneExp/'

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
		   
for cancer in cancers:
	normal = cancers[cancer]
	
	#compute human median TPM
	human_tpms = {}
	with open(tpm_dir+'human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms),2)
			human_tpms[gene] = [medianTpm]
	with open(tpm_dir+'human_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms),2)
			human_tpms[gene].append(medianTpm)
			
	#compute mouse median TPM
	mouse_tpms = {}
	with open(tpm_dir+'mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms),2)
			mouse_tpms[gene] = [medianTpm]
	with open(tpm_dir+'mouse_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			tpms = [float(e) for e in cols[1:]]
			medianTpm = round(np.median(tpms),2)
			mouse_tpms[gene].append(medianTpm)
	
	#read human normalized expression
	human_exps = {}
	with open(nx_dir+'human_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			human_exps[gene] = [exps]
			
	with open(nx_dir+'human_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			human_exps[gene].append(exps)
			
	#read mouse normalized expression
	mouse_exps = {}
	with open(nx_dir+'mouse_'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			mouse_exps[gene] = [exps]
			
	with open(nx_dir+'mouse_'+normal) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			exps = [float(e) for e in cols[1:]]
			mouse_exps[gene].append(exps)
			
	#two-sided Mann–Whitney test
	human_degs = []
	for gene in human_exps:
		if human_tpms[gene][0]>0 or human_tpms[gene][1]>0:
			medianExp1 = np.median(human_exps[gene][0])
			medianExp2 = np.median(human_exps[gene][1])
			deltaExp = round(medianExp1-medianExp2, 2)
			p = stats.mannwhitneyu(human_exps[gene][0], human_exps[gene][1])[1]
			human_degs.append([p, deltaExp, gene])
			
		
	mouse_degs = []
	for gene in mouse_exps:
		if mouse_tpms[gene][0]>0 or mouse_tpms[gene][1]>0:
			medianExp1 = np.median(mouse_exps[gene][0])
			medianExp2 = np.median(mouse_exps[gene][1])
			deltaExp = round(medianExp1-medianExp2, 2)
			p = stats.mannwhitneyu(mouse_exps[gene][0], mouse_exps[gene][1])[1]
			mouse_degs.append([p, deltaExp, gene])
		
	#Benjamini-Hochberg p-value correction
	human_pValues = [a[0] for a in human_degs]
	human_pValues = np.asfarray(human_pValues)
	by_descend = human_pValues.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(human_pValues)) / np.arange(len(human_pValues), 0, -1)
	fdrs = np.minimum(1, np.minimum.accumulate(steps * human_pValues[by_descend]))
	for i in range(len(fdrs[by_orig])):
		human_degs[i].append(fdrs[by_orig][i])
		
	mouse_pValues = [a[0] for a in mouse_degs]
	mouse_pValues = np.asfarray(mouse_pValues)
	by_descend = mouse_pValues.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(mouse_pValues)) / np.arange(len(mouse_pValues), 0, -1)
	fdrs = np.minimum(1, np.minimum.accumulate(steps * mouse_pValues[by_descend]))
	for i in range(len(fdrs[by_orig])):
		mouse_degs[i].append(fdrs[by_orig][i])
		
	#output DEG
	with open('DEG/human_'+cancer, 'w') as f:
		f.write('Gene\tDeltaNX\tCancerMedianLogTPM\tNormalMedianLogTPM\tPValue\tFDR\n')
		for p, deltaExp, gene, fdr in human_degs:
			cancermedianTPM, normalmedianTPM = human_tpms[gene]
			f.write('\t'.join([gene, str(deltaExp), str(cancermedianTPM), str(normalmedianTPM), str(p), str(fdr)])+'\n')
	
	with open('DEG/mouse_'+cancer, 'w') as f:
		f.write('Gene\tDeltaNX\tCancerMedianTPM\tNormalMedianTPM\tPValue\tFDR\n')
		for p, deltaExp, gene, fdr in mouse_degs:
			cancermedianTPM, normalmedianTPM = mouse_tpms[gene]
			f.write('\t'.join([gene, str(deltaExp), str(cancermedianTPM), str(normalmedianTPM), str(p), str(fdr)])+'\n')