#!/usr/bin/env python

from scipy import stats
import numpy as np
from sklearn.linear_model import LinearRegression
import math
import re

#linear regression
model = LinearRegression(copy_X=True,
						 fit_intercept=True,
						 n_jobs=5,
						 normalize=False)

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
		   
#cell types
cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
		 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
		 'NK cell','Monocyte','Tregs']

#read Ensembl IDs of orthologs
orthologs = {}
with open('gene_annotation/one-to-one-orthologs') as f:
	for line in f:
		ensg, ensm = line.strip().split('\t')
		orthologs[ensg] = ensg+'-'+ensm
		orthologs[ensm] = ensg+'-'+ensm
		
#extract co-expression modules
modules = {}
with open('co_expression_network/coExpModule') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		moduleID = cols[0]
		cancer = moduleID.split('|')[0]
		cancer = cancer.lower().replace(' ','_')
		
		moduleGenes = []
		humanGenes = cols[2].split(',')
		for gene in humanGenes:
			ensg = gene.split('|')[0]
			ensg_ensm = orthologs[ensg]
			moduleGenes.append(ensg_ensm)
			
		if cancer not in modules:
			modules[cancer] = {}
		modules[cancer][moduleID] = moduleGenes			
		
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
	
corr_results = {}	
for cancer in cancers:
	#read human gene expression
	gene_exps = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		samples = f.readline().strip().split('\t')[1:]
		samNum = len(samples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in orthologs:
				ensg_ensm = orthologs[cols[0]]
				exps = [float(e) for e in cols[1:]]
				gene_exps[ensg_ensm] = {}
				for i in range(samNum):
					gene_exps[ensg_ensm][samples[i]] = exps[i]
				
	#read mouse gene expression
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		mouseSamples = f.readline().strip().split('\t')[1:]
		samNum = len(mouseSamples)
		samples.extend(mouseSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in orthologs:
				ensg_ensm = orthologs[cols[0]]
				if ensg_ensm in gene_exps:
					exps = [float(e) for e in cols[1:]]
					for i in range(samNum):
						gene_exps[ensg_ensm][mouseSamples[i]] = exps[i]
	
	CIBERSORT_result = {}
	quanTIseq_result = {}
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
		mouseSamples = lines[0].strip().split(',')[1:]
		CIBERSORT_result['B cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0] in ('B cell naive_CIBERSORT-ABS','B cell memory_CIBERSORT-ABS','B cell plasma_CIBERSORT-ABS')]))
		CIBERSORT_result['T cell CD4+'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0] in ('T cell CD4+ naive_CIBERSORT-ABS','T cell CD4+ memory resting_CIBERSORT-ABS','T cell CD4+ memory activated_CIBERSORT-ABS')]))
		CIBERSORT_result['T cell CD8+'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='T cell CD8+_CIBERSORT-ABS']))
		CIBERSORT_result['Neutrophil'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Neutrophil_CIBERSORT-ABS']))
		CIBERSORT_result['Macrophage M1'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Macrophage M1_CIBERSORT-ABS']))
		CIBERSORT_result['Macrophage M2'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Macrophage M2_CIBERSORT-ABS']))
		CIBERSORT_result['Myeloid dendritic cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0] in ('Myeloid dendritic cell resting_CIBERSORT-ABS','Myeloid dendritic cell activated_CIBERSORT-ABS')]))
		CIBERSORT_result['Monocyte'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Monocyte_CIBERSORT-ABS']))
		CIBERSORT_result['NK cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0] in ('NK cell resting_CIBERSORT-ABS','NK cell activated_CIBERSORT-ABS')]))
		CIBERSORT_result['Tregs'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_CIBERSORT-ABS']))
		
		quanTIseq_result['B cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='B cell_QUANTISEQ']))
		quanTIseq_result['T cell CD4+'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='T cell CD4+ (non-regulatory)_QUANTISEQ']))
		quanTIseq_result['T cell CD8+'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='T cell CD8+_QUANTISEQ']))
		quanTIseq_result['Neutrophil'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Neutrophil_QUANTISEQ']))
		quanTIseq_result['Macrophage M1'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Macrophage M1_QUANTISEQ']))
		quanTIseq_result['Macrophage M2'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Macrophage M2_QUANTISEQ']))
		quanTIseq_result['Myeloid dendritic cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Myeloid dendritic cell_QUANTISEQ']))
		quanTIseq_result['Monocyte'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='Monocyte_QUANTISEQ']))
		quanTIseq_result['NK cell'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='NK cell_QUANTISEQ']))
		quanTIseq_result['Tregs'].update(sumLines(mouseSamples, [line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_QUANTISEQ']))
		
	#compute multiple correlation coefficient 
	#between immune infiltration and module genes
	for cell in cells:
		CIBERSORT_scores = []
		for sample in samples:
			if sample in CIBERSORT_result[cell]:
				CIBERSORT_scores.append(CIBERSORT_result[cell][sample])
		for moduleID in modules[cancer]:
			moduleExps = []
			moduleGenes = modules[cancer][moduleID]
			for sample in samples:
				if sample in CIBERSORT_result[cell]:
					moduleExps.append([gene_exps[ensg_ensm][sample] for ensg_ensm in moduleGenes])
			moduleExps_array = np.array(moduleExps)
			model.fit(moduleExps_array, CIBERSORT_scores)
			correlation = math.sqrt(model.score(moduleExps_array, CIBERSORT_scores))
			if moduleID not in corr_results:
				corr_results[moduleID] = {}
			corr_results[moduleID][cell+'|CIBERSORT'] = correlation
					
		quanTIseq_scores = []
		for sample in samples:
			if sample in quanTIseq_result[cell]:
				quanTIseq_scores.append(quanTIseq_result[cell][sample])
		for moduleID in modules[cancer]:
			moduleExps = []
			moduleGenes = modules[cancer][moduleID]
			for sample in samples:
				if sample in quanTIseq_result[cell]:
					moduleExps.append([gene_exps[ensg_ensm][sample] for ensg_ensm in moduleGenes])
			moduleExps_array = np.array(moduleExps)
			model.fit(moduleExps_array, quanTIseq_scores)
			correlation = math.sqrt(model.score(moduleExps_array, quanTIseq_scores))
			if moduleID not in corr_results:
				corr_results[moduleID] = {}
			corr_results[moduleID][cell+'|quanTIseq'] = correlation
				
#compute mean correlation
mean_corrs = {}
for moduleID in corr_results:
	mean_corrs[moduleID] = {}
	for cell in cells:
		mean_corr = round(np.mean((corr_results[moduleID][cell+'|CIBERSORT'], 
								   corr_results[moduleID][cell+'|quanTIseq'])),
						  2)
		mean_corrs[moduleID][cell] = mean_corr
		
#infer immune divergence-responsible modules
#extract immune infiltration divergence
diver_cells = {}
with open('immuneCell_comparison/immuneInfiltration_compare_result') as f:
	for line in f:
		cancer, cell = line.split('\t')[:2]
		cancer = cancer.capitalize()
		if cancer not in diver_cells:
			diver_cells[cancer] = set()
		diver_cells[cancer].add(cell)

#extract modules with high correlation with immune cell types		
module_genes = {}
with open('co_expression_network/coExpModule') as f: 
	f.readline()
	for line in f:
		cols = line.split('\t')
		moduleID = cols[0]
		cancer = moduleID.split('|')[0]
		for cell in mean_corrs[moduleID]:
			if mean_corrs[moduleID][cell] > 0.5 and cell in diver_cells[cancer]:
				module_genes[moduleID] = set(re.sub('ENSG\d+\|', '', cols[2]).split(','))
				break
				
#extract all genes in gtf file
gtf_genes = set()
with open('gene_annotation/Homo_sapiens.GRCh38.101.gtf') as f:
	for line in f:
		if '\tgene\t' in line:
			symbol = re.search('gene_name "(.+?)";', line).group(1)
			gtf_genes.add(symbol)
			
#extract immune response genes
immune_genes = set()
with open('gene_annotation/immune-GO/immune_response') as f:
	for line in f:
		symbol = line.strip().split('\t')[0]
		if symbol in gtf_genes:
			immune_genes.add(symbol)
			
#number of immune response genes
immuneGeneNum = len(immune_genes)

#extract all GO genes as background	for comparison
all_go_genes = set()
with open('data/go/goa_human.gaf') as f:
	for line in f:
		if line[0] != '!':
			symbol = line.split('\t')[2]
			if symbol in gtf_genes:
				all_go_genes.add(symbol)
				
#total number of GO genes
total_goGenes = len(all_go_genes)

#enrichment analysis via hypergeometric overlap test
enrichment = []
p_values = []	
for moduleID in module_genes:
	this_module_genes = module_genes[moduleID]
	
	hitCount = len(this_module_genes &
				   immune_genes)
	module_allGo_geneNum = len(this_module_genes &
							   all_go_genes)
	pVal = stats.hypergeom.sf(hitCount-1,total_goGenes,immuneGeneNum,module_allGo_geneNum)
	if pVal<0.05:
		enrichment.append([moduleID, pVal])
		p_values.append(pVal)
		
#Benjamini-Hochberg p-value correction
p_values = np.asfarray(p_values)
by_descend = p_values.argsort()[::-1]
by_orig = by_descend.argsort()
steps = float(len(p_values)) / np.arange(len(p_values), 0, -1)
fdrs = np.minimum(1, np.minimum.accumulate(steps * p_values[by_descend]))
for i in range(len(enrichment)):
	enrichment[i].append(fdrs[by_orig][i])
	
#immune divergence-responsible modules
immune_modules = set()
for moduleID, pVal, fdr in enrichment:
	if fdr<0.05:
		immune_modules.add(moduleID)

#write results
f_re = open('co_expression_network/coExpModule_immuneInfiltration', 'w')
with open('co_expression_network/coExpModule') as f:
	title_cols = f.readline().split('\t')
	newTitle_cols = title_cols[:2]+\
					['Immune divergence-responsible module(Yes or No)']+\
					cells+\
					['Divergent genes(human symbol)']+\
					title_cols[4:]
	f_re.write('\t'.join(newTitle_cols))
	for line in f:
		cols = line.split('\t')
		if cols[0] in immune_modules:
			respo_label = 'Y'
		else:
			respo_label = 'N'
		newCols = cols[:2]+\
				  [respo_label]+\
				  [str(mean_corrs[cols[0]][cell]) for cell in cells]+\
				  [re.sub('ENSG\d+\|', '', cols[2])]+\
				  [cols[4]]+\
				  [re.sub('ENSG\d+\|', '', cols[5])]+\
				  [cols[6]]+\
				  [re.sub('ENSG\d+\|', '', cols[7])]+\
				  [cols[8]]+\
				  [re.sub('ENSMUSG\d+\|', '', cols[9])]+\
				  [cols[10]]
		f_re.write('\t'.join(newCols))
f_re.close()