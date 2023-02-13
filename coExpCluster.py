#!/usr/bin/env python

from sklearn.cluster import AgglomerativeClustering
from scipy import stats
import pandas as pd
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
		   
#extract differential expression results of human VS mouse
de_results = {}
for cancer in cancers:
	de_results[cancer] = {}
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			deltaNX = float(cols[4])
			fdr = float(cols[-1])
			if abs(deltaNX)>1 and fdr<0.05:
				de_results[cancer][cols[0]] = {}
				de_results[cancer][cols[1]] = {}
				if deltaNX>0:
					de_results[cancer][cols[0]]['human'] = 'up'
					de_results[cancer][cols[1]]['mouse'] = 'down'
				else:
					de_results[cancer][cols[0]]['human'] = 'down'
					de_results[cancer][cols[1]]['mouse'] = 'up'
		   
#extract lncRNAs' clade-specificity
specificity = {}
geneSymbols = {}
with open('lncrna_analysis/primateSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		specificity[cols[0]] = cols[6]
		geneSymbols[cols[0]] = cols[0]+'|'+cols[1]
with open('lncrna_analysis/rodentSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		specificity[cols[0]] = cols[6]
		geneSymbols[cols[0]] = cols[0]+'|'+cols[1]
		
#extract symbols of divergent genes
for cancer in cancers:
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			ensg, ensm, symbolH, symbolM = line.split('\t')[:4]
			geneSymbols[ensg] = ensg+'|'+symbolH
			geneSymbols[ensm] = ensm+'|'+symbolM

f_re = open('co_expression_network/coExpModule', 'w')
f_re2 = open('co_expression_network/coExpModule_lncExample', 'w')
f_re.write('\t'.join(['Module',
					  'Size',
					  'DiverGenes(HumanID)',
					  'DiverGenes(MouseID)',
					  'DiverGene_num',
					  'SimianLnc',
					  'SimianLnc_num',
					  'PrimateLnc',
					  'PrimateLnc_num',
					  'RodentLnc',
					  'RodentLnc_num'])+'\n')
f_re2.write('\t'.join(['LncRNA',
					   'Module',
					   'FDR',
					   'MedianAbscorrModule',
					   'MedianAbscorrOther',
					   'ModuleGeneNum',
					   'OtherGeneNum'])+'\n')					  
for cancer in cancers:
	#compute correlation distance
	distances = {}
	with open('co_expression_network/divergentGene_correlation/'+cancer) as f:
		for line in f:
			gene1, gene2, corr, p = line.strip().split('\t')
			if gene1 not in distances:
				distances[gene1] = {}
				distances[gene1][gene1] = 0
			if gene2 not in distances:
				distances[gene2] = {}
				distances[gene2][gene2] = 0
			if gene2 not in distances[gene1]:
				distances[gene1][gene2] = 1-abs(float(corr))
			if gene1 not in distances[gene2]:
				distances[gene2][gene1] = 1-abs(float(corr))
		
	#convert to dataframe
	distances_df = pd.DataFrame(distances)
	
	#distance array
	distance_array = distances_df.values
	
	#agglomerative clustering
	model = AgglomerativeClustering(n_clusters=None,
									affinity='precomputed',
									linkage='complete',
									distance_threshold=np.percentile(distance_array,75),
									compute_full_tree=True
									)
	clustering = model.fit(distance_array)
	clusterLabels = list(clustering.labels_)
	
	#generate modules
	pre_modules = {}
	diverGenes = list(distances_df.index)
	for i in range(len(diverGenes)):
		moduleNum = clusterLabels[i]
		moduleName = cancer.replace('_',' ').capitalize()+'|ME'+str(moduleNum)
		diverGene_h, diverGene_m = diverGenes[i].split('|')
		if moduleName not in pre_modules:
			pre_modules[moduleName] = set()
		pre_modules[moduleName].add(diverGene_h)
		pre_modules[moduleName].add(diverGene_m)
		
	#abandon modules with size<=5
	modules = {}
	num = 0
	for moduleName in pre_modules:
		size = len([gene for gene in pre_modules[moduleName] if gene[:4]=='ENSG'])
		if size >5 :
			new_moduleName = cancer.replace('_',' ').capitalize()+'|ME'+str(num)
			modules[new_moduleName] = pre_modules[moduleName]
			num += 1
	
	#extract lncRNA-divergent gene correlation
	correlation = {}
	with open('co_expression_network/cladeLnc_divergentGene_correlation/human_'+cancer) as f:
		for line in f:
			lncID, diverGeneID, corr, p = line.split('\t')
			if lncID not in correlation:
				correlation[lncID] = {}
			if de_results[cancer][diverGeneID]['human'] == 'up':
				correlation[lncID][diverGeneID] = float(corr)
			elif de_results[cancer][diverGeneID]['human'] == 'down':
				correlation[lncID][diverGeneID] = -float(corr)
	with open('co_expression_network/cladeLnc_divergentGene_correlation/mouse_'+cancer) as f:
		for line in f:
			lncID, diverGeneID, corr, p = line.split('\t')
			if lncID not in correlation:
				correlation[lncID] = {}
			if de_results[cancer][diverGeneID]['mouse'] == 'up':
				correlation[lncID][diverGeneID] = float(corr)
			elif de_results[cancer][diverGeneID]['mouse'] == 'down':
				correlation[lncID][diverGeneID] = -float(corr)

	#add lncRNA to module using two-sample Kolmogorov-Smirnov test
	module_lncs = {}
	module_lncs_list = []
	pValues = []
	lncMedianCorrs = {}
	for moduleName in modules:
		for lncID in correlation:
			lncModule_corrs = [correlation[lncID][diverGeneID] for 
							   diverGeneID in modules[moduleName] 
							   if diverGeneID in correlation[lncID]]
			lncOther_corrs = [correlation[lncID][geneID] for geneID 
							   in correlation[lncID] if geneID not 
							   in modules[moduleName]]
			p = stats.ks_2samp(lncModule_corrs, lncOther_corrs, alternative='less')[1]
			module_lncs_list.append((lncID, moduleName))
			pValues.append(p)
			lncMedianCorrs[(lncID, moduleName)] = [np.median(lncModule_corrs),
												   np.median(lncOther_corrs),
												   len(lncModule_corrs),
												   len(lncOther_corrs)]
			
	#FDR correction
	pValues = np.asfarray(pValues)
	by_descend = pValues.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
	fdrs = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
	fdrNum = len(fdrs[by_orig])
	for i in range(fdrNum):
		if fdrs[by_orig][i] < 0.00001:
			lncID, moduleName = module_lncs_list[i]
			if moduleName not in module_lncs:
				module_lncs[moduleName] = set()
			module_lncs[moduleName].add(lncID)
			
			lncSymbols = geneSymbols[lncID]
			median1, median2, num1, num2 = lncMedianCorrs[(lncID, moduleName)]
			f_re2.write('\t'.join([lncSymbols,
								   moduleName,
								   str(fdrs[by_orig][i]),
								   str(median1),
								   str(median2),
								   str(num1),
								   str(num2)])+'\n')
			
	for moduleName in modules:
		genes_h = [geneSymbols[geneID] for geneID in modules[moduleName] if geneID[:4]=='ENSG']
		genes_m = [geneSymbols[geneID] for geneID in modules[moduleName] if geneID[:4]=='ENSM']
		diverNum = len(genes_h)
		simianLncs = []
		primateLncs = []
		rodentLncs = []
		if moduleName in module_lncs:
			for lncID in module_lncs[moduleName]:
				lnc_speci = specificity[lncID]
				if lnc_speci in ('simian-specific','human-specific'):
					simianLncs.append(geneSymbols[lncID])
				elif lnc_speci == 'primate-specific':
					primateLncs.append(geneSymbols[lncID])
				elif lnc_speci == 'rodent-specific':
					rodentLncs.append(geneSymbols[lncID])
					
		simianNum = len(simianLncs)
		primateNum = len(primateLncs)
		rodentNum = len(rodentLncs)
		size = diverNum+simianNum+primateNum+rodentNum
		
		genes_h_str = ','.join(genes_h)
		genes_m_str = ','.join(genes_m)
		simianLncs_str = ','.join(simianLncs)
		primateLncs_str = ','.join(primateLncs)
		rodentLncs_str = ','.join(rodentLncs)
		f_re.write('\t'.join([moduleName,
						   str(size),
						   genes_h_str,
						   genes_m_str,
						   str(diverNum),
						   simianLncs_str,
						   str(simianNum),
						   primateLncs_str,
						   str(primateNum),
						   rodentLncs_str,
						   str(rodentNum)])+'\n')
f_re.close()
f_re2.close()