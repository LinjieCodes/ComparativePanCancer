#!/usr/bin/env python

from scipy import stats
import numpy as np
import os
import re

allGenes = set()
with open('gene_annotation/Homo_sapiens.GRCh38.101.gtf') as f:
	for line in f:
		if '\tgene\t' in line and 'gene_biotype "protein_coding";' in line:
			gene = re.search('gene_name "(.+?)";', line).group(1)
			allGenes.add(gene)
allGeneNum = len(allGenes)

markers = {}
markerFilePath = 'gene_annotation/CellMarker/'
for file in os.listdir(markerFilePath):
	cell = file.split('_')[1]
	this_markers = set()
	with open(markerFilePath+file) as f:
		for line in f:
			cols = line.split('\t')
			genes = cols[4].split(', ')
			for gene in genes:
				geneUpp = gene.upper()
				if geneUpp in allGenes:
					this_markers.add(geneUpp)
	if this_markers:
		markers[cell] = this_markers
	else:
		print('Markers of ', cell, ': None!')

enrichment = []
p_values = []
with open('coExpModule_immuneInfiltration') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		if cols[2] == 'Y':
			module = cols[0]
			module_genes = set(cols[13].split(','))
			module_gene_num = len(module_genes)
			for cell in markers:
				hitCount = len(module_genes &
							   markers[cell])
				this_marker_num = len(markers[cell])
				pVal = stats.hypergeom.sf(hitCount - 1,
										  allGeneNum,
										  module_gene_num,
										  this_marker_num)
				if pVal<0.05:
					enrichment.append([module, cell, pVal])
					p_values.append(pVal)
					
#Benjamini-Hochberg p-value correction
p_values = np.asfarray(p_values)
by_descend = p_values.argsort()[::-1]
by_orig = by_descend.argsort()
steps = float(len(p_values)) / np.arange(len(p_values), 0, -1)
fdrs = np.minimum(1, np.minimum.accumulate(steps * p_values[by_descend]))
for i in range(len(enrichment)):
	enrichment[i].append(fdrs[by_orig][i])
	
with open('MarkerEnrichment', 'w') as f:
	for module, cell, p, fdr in enrichment:
		if fdr<0.1:
			f.write('\t'.join([module,
							   cell,
							   str(p),
							   str(fdr)])+'\n')