#!/usr/bin/env python

from scipy import stats
import numpy as np

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
 
#get significantly divergent cell types
score_fc = {}
cancer_cells = {}
for cancer in cancers:
	cancer_cells[cancer] = {}
with open('immuneCell_comparison/immuneInfiltration_compare_result') as f:
	for line in f:
		cols = line.split('\t')
		cancer_cells[cols[0].replace(' ','_')][cols[1]] = []
		cancer_cells[cols[0].replace(' ','_')][cols[1]].append('CIBERSORT')
		cancer_cells[cols[0].replace(' ','_')][cols[1]].append('quanTIseq')
			
#read module genes
module_lncs_h = {}
module_lncs_m = {}
module_genes = {}
with open('co_expression_network/coExpModule_immuneInfiltration') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		if cols[2] == 'Y':
			module = cols[0]
			genes = set(cols[13].split(','))
			humanLncs = set(cols[15].split(',')) | set(cols[17].split(','))
			mouseLncs = set(cols[19].split(','))
			if humanLncs != set(['']):
				if '' in humanLncs:
					humanLncs.remove('')
				module_lncs_h[module] = humanLncs
			if mouseLncs != set(['']):
				module_lncs_m[module] = mouseLncs
			module_genes[module] = genes

def sumLines(mouseSamples, lines):
	sumList = [0]*len(mouseSamples)
	for line in lines:
		thisList = [float(value) for value in line.strip().split(',')[1:]]
		for index in range(len(thisList)):
			sumList[index]+= thisList[index]
	this_infiltration = {}
	for i in range(len(mouseSamples)):
		this_infiltration[mouseSamples[i]] = sumList[i]
	return this_infiltration
	
#get lncRNA ENSG-symbol
geneSymbols = {}
with open('lncrna_analysis/primateSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		geneSymbols[cols[0]] = cols[1]
with open('lncrna_analysis/rodentSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		geneSymbols[cols[0]] = cols[1]

results_cbs = {}
results_qts = {}
p_values_cbs = []
p_values_qts = []
geneCancerCells_cbs = []
geneCancerCells_qts = []	
for cancer in cancers:
	#get divergent genes ENSG-symbol
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			if abs(float(cols[4]))>1 and float(cols[-1])<0.05:
				geneSymbols[cols[0]] = cols[2]
				geneSymbols[cols[1]] = cols[2]
		
	#read human gene expression
	gene_exps = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		humanSamples = f.readline().strip().split('\t')[1:]
		samNum = len(humanSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneSymbols:
				symbol = geneSymbols[cols[0]]
				exps = [float(e) for e in cols[1:]]
				gene_exps[symbol] = {}
				for i in range(samNum):
					gene_exps[symbol][humanSamples[i]] = exps[i]
				
	#read mouse gene expression
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		mouseSamples = f.readline().strip().split('\t')[1:]
		samNum = len(mouseSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneSymbols:
				symbol = geneSymbols[cols[0]]
				exps = [float(e) for e in cols[1:]]
				if symbol not in gene_exps:	#indicates mouse lncRNA
					gene_exps[symbol] = {}
				for i in range(samNum):
					gene_exps[symbol][mouseSamples[i]] = exps[i]
					
	#all samples
	samples = humanSamples + mouseSamples
	
	CIBERSORT_result = {}
	quanTIseq_result = {}
	cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
			 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
			 'NK cell','Monocyte','Tregs']
	for cell in cells:
		CIBERSORT_result[cell] = {}
		quanTIseq_result[cell] = {}
	#read infiltration estimation of human
	with open('immuneCell_comparison/TIMER2.0/infiltration_estimation_for_tcga.csv') as f:
		for line in f:
			cols = line.strip().split(',')
			thisSample = cols[0]
			if thisSample in samples:
				values = [float(value) for value in cols[1:] if value !='NA']
				if cancer=='leukemia':#TIMER score is not available for leukemia
					values = ['NA','NA','NA','NA','NA','NA']+values
				CIBERSORT_result['B cell'][thisSample] = sum(values[28:31])
				CIBERSORT_result['T cell CD4+'][thisSample] = sum(values[32:35])
				CIBERSORT_result['T cell CD8+'][thisSample] = values[31]
				CIBERSORT_result['Neutrophil'][thisSample] = values[49]
				CIBERSORT_result['Macrophage M1'][thisSample] = values[42]
				CIBERSORT_result['Macrophage M2'][thisSample] = values[43]
				CIBERSORT_result['Myeloid dendritic cell'][thisSample] = sum(values[44:46])
				CIBERSORT_result['Monocyte'][thisSample] = values[40]
				CIBERSORT_result['NK cell'][thisSample] = sum(values[38:40])
				CIBERSORT_result['Tregs'][thisSample] = values[36]
				
				quanTIseq_result['B cell'][thisSample] = values[50]
				quanTIseq_result['T cell CD4+'][thisSample] = values[56]
				quanTIseq_result['T cell CD8+'][thisSample] = values[57]
				quanTIseq_result['Neutrophil'][thisSample] = values[54]
				quanTIseq_result['Macrophage M1'][thisSample] = values[51]
				quanTIseq_result['Macrophage M2'][thisSample] = values[52]
				quanTIseq_result['Myeloid dendritic cell'][thisSample] = values[59]
				quanTIseq_result['Monocyte'][thisSample] = values[53]
				quanTIseq_result['NK cell'][thisSample] = values[55]
				quanTIseq_result['Tregs'][thisSample] = values[58]
			
	#read infiltration estimation of mouse
	with open('immuneCell_comparison/TIMER2.0/mouse_%s.csv' % cancer) as f:
		lines = f.readlines()
		mouseSamples_cell = lines[0].strip().split(',')[1:]
		CIBERSORT_result['B cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0] in ('B cell naive_CIBERSORT-ABS','B cell memory_CIBERSORT-ABS','B cell plasma_CIBERSORT-ABS')]))
		CIBERSORT_result['T cell CD4+'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0] in ('T cell CD4+ naive_CIBERSORT-ABS','T cell CD4+ memory resting_CIBERSORT-ABS','T cell CD4+ memory activated_CIBERSORT-ABS')]))
		CIBERSORT_result['T cell CD8+'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='T cell CD8+_CIBERSORT-ABS']))
		CIBERSORT_result['Neutrophil'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Neutrophil_CIBERSORT-ABS']))
		CIBERSORT_result['Macrophage M1'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Macrophage M1_CIBERSORT-ABS']))
		CIBERSORT_result['Macrophage M2'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Macrophage M2_CIBERSORT-ABS']))
		CIBERSORT_result['Myeloid dendritic cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0] in ('Myeloid dendritic cell resting_CIBERSORT-ABS','Myeloid dendritic cell activated_CIBERSORT-ABS')]))
		CIBERSORT_result['Monocyte'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Monocyte_CIBERSORT-ABS']))
		CIBERSORT_result['NK cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0] in ('NK cell resting_CIBERSORT-ABS','NK cell activated_CIBERSORT-ABS')]))
		CIBERSORT_result['Tregs'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_CIBERSORT-ABS']))
		
		quanTIseq_result['B cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='B cell_QUANTISEQ']))
		quanTIseq_result['T cell CD4+'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='T cell CD4+ (non-regulatory)_QUANTISEQ']))
		quanTIseq_result['T cell CD8+'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='T cell CD8+_QUANTISEQ']))
		quanTIseq_result['Neutrophil'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Neutrophil_QUANTISEQ']))
		quanTIseq_result['Macrophage M1'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Macrophage M1_QUANTISEQ']))
		quanTIseq_result['Macrophage M2'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Macrophage M2_QUANTISEQ']))
		quanTIseq_result['Myeloid dendritic cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Myeloid dendritic cell_QUANTISEQ']))
		quanTIseq_result['Monocyte'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='Monocyte_QUANTISEQ']))
		quanTIseq_result['NK cell'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='NK cell_QUANTISEQ']))
		quanTIseq_result['Tregs'].update(sumLines(mouseSamples_cell, [line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_QUANTISEQ']))
		
	#compute Spearman correlation between infiltration and gene expression
	for cell in cancer_cells[cancer]:
		if 'CIBERSORT' in cancer_cells[cancer][cell]:
			#divergent gene, CIBERSORT
			CIBERSORT_scores = []
			for sample in samples:
				if sample in CIBERSORT_result[cell]:
					CIBERSORT_scores.append(CIBERSORT_result[cell][sample])
			for module in module_genes:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_genes[module]:
					thisGeneExps = []
					for sample in samples:
						if sample in CIBERSORT_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(CIBERSORT_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_cbs.append(pVal)
					geneCancerCells_cbs.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_cbs:
						results_cbs[(symbol, module)] = {}
					results_cbs[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
					
			#human lncRNAs, CIBERSORT
			CIBERSORT_scores = []
			for sample in humanSamples:
				if sample in CIBERSORT_result[cell]:
					CIBERSORT_scores.append(CIBERSORT_result[cell][sample])
			for module in module_lncs_h:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_lncs_h[module]:
					thisGeneExps = []
					for sample in humanSamples:
						if sample in CIBERSORT_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(CIBERSORT_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_cbs.append(pVal)
					geneCancerCells_cbs.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_cbs:
						results_cbs[(symbol, module)] = {}
					results_cbs[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
					
			#mouse lncRNAs, CIBERSORT
			CIBERSORT_scores = []
			for sample in mouseSamples:
				if sample in CIBERSORT_result[cell]:
					CIBERSORT_scores.append(CIBERSORT_result[cell][sample])
			for module in module_lncs_m:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_lncs_m[module]:
					thisGeneExps = []
					for sample in mouseSamples:
						if sample in CIBERSORT_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(CIBERSORT_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_cbs.append(pVal)
					geneCancerCells_cbs.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_cbs:
						results_cbs[(symbol, module)] = {}
					results_cbs[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
					
		if 'quanTIseq' in cancer_cells[cancer][cell]:
			#divergent gene, quanTIseq
			quanTIseq_scores = []
			for sample in samples:
				if sample in quanTIseq_result[cell]:
					quanTIseq_scores.append(quanTIseq_result[cell][sample])
			for module in module_genes:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_genes[module]:
					thisGeneExps = []
					for sample in samples:
						if sample in quanTIseq_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(quanTIseq_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_qts.append(pVal)
					geneCancerCells_qts.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_qts:
						results_qts[(symbol, module)] = {}
					results_qts[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
					
			#human lncRNAs, quanTIseq
			quanTIseq_scores = []
			for sample in humanSamples:
				if sample in quanTIseq_result[cell]:
					quanTIseq_scores.append(quanTIseq_result[cell][sample])
			for module in module_lncs_h:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_lncs_h[module]:
					thisGeneExps = []
					for sample in humanSamples:
						if sample in quanTIseq_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(quanTIseq_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_qts.append(pVal)
					geneCancerCells_qts.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_qts:
						results_qts[(symbol, module)] = {}
					results_qts[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
					
			#mouse gene, quanTIseq
			quanTIseq_scores = []
			for sample in mouseSamples:
				if sample in quanTIseq_result[cell]:
					quanTIseq_scores.append(quanTIseq_result[cell][sample])
			for module in module_lncs_m:
				module_cancer = module.split('|')[0].replace(' ','_').lower()
				if module_cancer != cancer:
					continue
				for symbol in module_lncs_m[module]:
					thisGeneExps = []
					for sample in mouseSamples:
						if sample in quanTIseq_result[cell]:
							thisGeneExps.append(gene_exps[symbol][sample])
					correlation,pVal = stats.spearmanr(quanTIseq_scores, thisGeneExps)
					if np.isnan(pVal):
						continue
					p_values_qts.append(pVal)
					geneCancerCells_qts.append((symbol, module, cancer.replace('_',' ')+'|'+cell))
					if (symbol, module) not in results_qts:
						results_qts[(symbol, module)] = {}
					results_qts[(symbol, module)][cancer.replace('_',' ')+'|'+cell] = round(correlation,2)
				
#Benjamini-Hochberg p-value correction
significant_correlations_cbs = set()
p_values_cbs = np.asfarray(p_values_cbs)
by_descend = p_values_cbs.argsort()[::-1]
by_orig = by_descend.argsort()
steps = float(len(p_values_cbs)) / np.arange(len(p_values_cbs), 0, -1)
fdrs = np.minimum(1, np.minimum.accumulate(steps * p_values_cbs[by_descend]))
for i in range(len(geneCancerCells_cbs)):
	if fdrs[by_orig][i]<0.05:
		significant_correlations_cbs.add((geneCancerCells_cbs[i], fdrs[by_orig][i]))
		
significant_correlations_qts = set()
p_values_qts = np.asfarray(p_values_qts)
by_descend = p_values_qts.argsort()[::-1]
by_orig = by_descend.argsort()
steps = float(len(p_values_qts)) / np.arange(len(p_values_qts), 0, -1)
fdrs = np.minimum(1, np.minimum.accumulate(steps * p_values_qts[by_descend]))
for i in range(len(geneCancerCells_qts)):
	if fdrs[by_orig][i]<0.05:
		significant_correlations_qts.add((geneCancerCells_qts[i], fdrs[by_orig][i]))

#write results
with open('Spearman_CIBERSORTx', 'w') as f:
	f.write('\t'.join(['Symbol',
					   'Module',
					   'Cell',
					   'rho',
					   'FDR'])+'\n')
	for t, fdr in significant_correlations_cbs:
		symbol, module, cell = t
		rho = results_cbs[(symbol, module)][cell]
		f.write('\t'.join([symbol,
						   module,
						   cell,
						   str(rho),
						   str(fdr)])+'\n')
						   
with open('Spearman_quanTIseq', 'w') as f:
	f.write('\t'.join(['Symbol',
					   'Module',
					   'Cell',
					   'rho',
					   'FDR'])+'\n')
	for t, fdr in significant_correlations_qts:
		symbol, module, cell = t
		rho = results_qts[(symbol, module)][cell]
		f.write('\t'.join([symbol,
						   module,
						   cell,
						   str(rho),
						   str(fdr)])+'\n')