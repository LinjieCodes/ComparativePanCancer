#!/usr/bin/env python

import math

#extract FDRs of KS test
fdrs = {}
lncrnas = set()
with open('co_expression_network/coExpModule_lncExample') as f:
	f.readline()
	for line in f:
		lncID_sym, module, fdr = line.split('\t')[:3]
		fdrs[(lncID_sym, module)] = -math.log10(float(fdr))
		lncrnas.add(lncID_sym)
		
#extract immune divergence-responsible modules
immune_modules = {}
with open('co_expression_network/coExpModule_immuneInfiltration') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		if cols[2] == 'Y':
			r_list = [float(r) for r in cols[3:12]]
			max_r = max(r_list)
			immune_modules[cols[0]] = max_r

#13 cancer types
cancers = ['Bladder cancer',
		   'Breast cancer',
		   'Glioblastoma',
		   'Glioma',
		   'Kidney cancer',
		   'Leukemia',
		   'Liver cancer',
		   'Lung cancer',
		   'Lymphoma',
		   'Ovarian cancer',
		   'Pancreatic cancer',
		   'Prostate cancer',
		   'Skin melanoma']		   
	
#calculate IR scores
ir_scores = {}
lncrna_modules = {}
for lncID_sym in lncrnas:
	ir_scores[lncID_sym] = {}
	lncrna_modules[lncID_sym] = set()
	for module in immune_modules:
		if (lncID_sym, module) in fdrs:
			fdr = fdrs[(lncID_sym, module)]
			max_r = immune_modules[module]
			ir_socre = fdr * max_r
			cancer = module.split('|')[0]
			if cancer not in ir_scores[lncID_sym]:
				ir_scores[lncID_sym][cancer] = []
			ir_scores[lncID_sym][cancer].append(ir_socre)
			lncrna_modules[lncID_sym].add(module)
			
#calculate sum IR scores
sir_scores_human = []
sir_scores_mouse = []
cancerNums = {}
for lncID_sym in ir_scores:
	sir_score = 0
	cancerNums[lncID_sym] = 0
	for cancer in cancers:
		if cancer in ir_scores[lncID_sym]:
			sir_score += max(ir_scores[lncID_sym][cancer])
			cancerNums[lncID_sym] += 1
		else:
			sir_score += 0
	if lncID_sym[:4]=='ENSG':
		sir_scores_human.append((sir_score, lncID_sym))
	else:
		sir_scores_mouse.append((sir_score, lncID_sym))
		
#read lncRNAs' lineage-specificity
lncrna_lineage = {}
with open('co_expression_network/coExpModule_immuneInfiltration') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		if cols[15]:
			simianLncs = cols[15].split(',')
			for lnc_sym in simianLncs:
				lncrna_lineage[lnc_sym] = 'simian-specific'
		if cols[17]:
			primateLncs = cols[17].split(',')
			for lnc_sym in primateLncs:
				lncrna_lineage[lnc_sym] = 'primate-specific'
		if cols[19]:
			rodentLncs = cols[19].split(',')
			for lnc_sym in rodentLncs:
				lncrna_lineage[lnc_sym] = 'rodent-specific'

#top10% of lncRNAs				
sir_scores_human.sort(reverse=True)
sir_scores_mouse.sort(reverse=True)
with open('co_expression_network/immune_divergence_responsible_lncRNAs','w') as f:
	f.write('\t'.join(['LncRNA symbol',
					   'Ensembl ID',
					   'Species',
					   'Lineage-specificity',
					   'Sum of IR scores',
					   'Number of cancers in which the lncRNA regulates immune divergence-responsible modules',
					   'In which immune divergence-responsible modules'])+'\n')

	for score, lncID_sym in sir_scores_human[:int(len(sir_scores_human)*0.1)]:
		lncID, lnc_sym = lncID_sym.split('|')
		species = 'Human'
		if lnc_sym not in lncrna_lineage:
			continue
		specificity = lncrna_lineage[lnc_sym]
		cancerNum = cancerNums[lncID_sym]
		modules_str = ','.join(lncrna_modules[lncID_sym])
		f.write('\t'.join([lnc_sym,
						   lncID,
						   species,
						   specificity,
						   str(round(score,2)),
						   str(cancerNum),
						   modules_str])+'\n')

	for score, lncID_sym in sir_scores_mouse[:int(len(sir_scores_mouse)*0.1)]:
		lncID, lnc_sym = lncID_sym.split('|')
		species = 'Mouse'
		if lnc_sym not in lncrna_lineage:
			continue
		specificity = lncrna_lineage[lnc_sym]
		cancerNum = cancerNums[lncID_sym]
		modules_str = ','.join(lncrna_modules[lncID_sym])
		f.write('\t'.join([lnc_sym,
						   lncID,
						   species,
						   specificity,
						   str(round(score,2)),
						   str(cancerNum),
						   modules_str])+'\n')