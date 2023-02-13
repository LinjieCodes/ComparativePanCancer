#!/usr/bin/env python

import re
import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams['font.sans-serif']=['Arial','sans-serif']
		
#read primate- and rodent-specific lncRNAs
cladeLncs_primate = set()
cladeLncs_rodent = set()
conservations = {}
with open('gene_annotation/primateSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		ensg = cols[0]
		conser = cols[6]
		cladeLncs_primate.add(ensg)
		conservations[ensg] = conser
with open('gene_annotation/rodentSpecific_lncRNA_summary') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		ensm = cols[0]
		conser = cols[6]
		cladeLncs_rodent.add(ensm)
		conservations[ensm] = conser
		
#read conserved lncRNAs
conserLncs_human = set()
conserLncs_mouse = set()
with open('gene_annotation/gencode.v36.long_noncoding_RNAs.gtf') as f:
	for line in f:
		if '\tgene\t' in line:
			ensg = re.search('gene_id "(.+?)"', line).group(1)
			ensg = ensg.split('.')[0]
			if ensg not in cladeLncs_primate:
				conserLncs_human.add(ensg)
with open('gene_annotation/gencode.vM22.long_noncoding_RNAs.gtf') as f:
	for line in f:
		if '\tgene\t' in line:
			ensm = re.search('gene_id "(.+?)"', line).group(1)
			ensm = ensm.split('.')[0]
			if ensm not in cladeLncs_rodent:
				conserLncs_mouse.add(ensm)
				
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
#retain lncRNAs and TFs with TPM>0.1 in 
#>=50% of samples in each cancer group
conserLncs_human_expressed = {}
conserLncs_mouse_expressed = {}
cladeLncs_primate_expressed = {}
cladeLncs_rodent_expressed = {}
medianTPM = {}
for cancer in cancers:
	medianTPM[cancer] = {}
	with open('data/geneTPM_ComBatCorrected_log2/human_'+cancer) as f:
		sampleNum = len(f.readline().split('\t'))-1
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in conserLncs_human:
				exps = [2**float(e) for e in cols[1:]]
				expressed_sampleNum = len([e for e in exps if e>0.1])
				if expressed_sampleNum>=0.5*sampleNum:
					if cancer not in conserLncs_human_expressed:
						conserLncs_human_expressed[cancer] = set()
					conserLncs_human_expressed[cancer].add(cols[0])
					medianTPM[cancer][cols[0]] = round(np.median(exps),2)
			elif cols[0] in cladeLncs_primate:
				exps = [2**float(e) for e in cols[1:]]
				expressed_sampleNum = len([e for e in exps if e>0.1])
				if expressed_sampleNum>=0.5*sampleNum:
					if cancer not in cladeLncs_primate_expressed:
						cladeLncs_primate_expressed[cancer] = set()
					cladeLncs_primate_expressed[cancer].add(cols[0])
					medianTPM[cancer][cols[0]] = round(np.median(exps),2)
	with open('data/geneTPM_ComBatCorrected_log2/mouse_'+cancer) as f:
		sampleNum = len(f.readline().split('\t'))-1
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in conserLncs_mouse:
				exps = [2**float(e) for e in cols[1:]]
				expressed_sampleNum = len([e for e in exps if e>0.1])
				if expressed_sampleNum>=0.5*sampleNum:
					if cancer not in conserLncs_mouse_expressed:
						conserLncs_mouse_expressed[cancer] = set()
					conserLncs_mouse_expressed[cancer].add(cols[0])
					medianTPM[cancer][cols[0]] = round(np.median(exps),2)
			elif cols[0] in cladeLncs_rodent:
				exps = [2**float(e) for e in cols[1:]]
				expressed_sampleNum = len([e for e in exps if e>0.1])
				if expressed_sampleNum>=0.5*sampleNum:
					if cancer not in cladeLncs_rodent_expressed:
						cladeLncs_rodent_expressed[cancer] = set()
					cladeLncs_rodent_expressed[cancer].add(cols[0])
					medianTPM[cancer][cols[0]] = round(np.median(exps),2)
				
#write result
numbers = []
for cancer in cancers:
	simianLncNum = len([lnc for lnc in cladeLncs_primate_expressed[cancer] if 
	                    conservations[lnc]=='simian-specific'])
	primateLncNum = len(cladeLncs_primate_expressed[cancer]) - simianLncNum
	humanConserLncNum = len(conserLncs_human_expressed[cancer])
	humanTotal = humanConserLncNum+simianLncNum+primateLncNum
	
	rodentLncNum = len(cladeLncs_rodent_expressed[cancer])
	mouseConcerLncNum = len(conserLncs_mouse_expressed[cancer])
	mouseTotal = mouseConcerLncNum+rodentLncNum
	print('##############')
	print(cancer)
	print('simian-specific:', simianLncNum)
	print('primate-specific:', primateLncNum)
	print('Total human lncRNAs:', humanTotal)
	print()
	print('rodent-specific:', rodentLncNum)
	print('Total mouse lncRNAs:', mouseTotal)
	print()
	numbers.append((cancer, simianLncNum, primateLncNum, humanConserLncNum, rodentLncNum, mouseConcerLncNum))
	
	
	with open('lncrna_analysis/cladeLnc_inCancers/human_'+cancer, 'w') as f:
		f.write('LncRNA\tClade-specificity\tMedianTPM\n')
		for lnc in cladeLncs_primate_expressed[cancer]:
			f.write('\t'.join([lnc, conservations[lnc], str(medianTPM[cancer][lnc])])+'\n')
			
	with open('lncrna_analysis/cladeLnc_inCancers/mouse_'+cancer, 'w') as f:
		f.write('LncRNA\tClade-specificity\tMedianTPM\n')
		for lnc in cladeLncs_rodent_expressed[cancer]:
			f.write('\t'.join([lnc, conservations[lnc], str(medianTPM[cancer][lnc])])+'\n')
			
#plot numbers
width=0.1
human_xs = [0.5+2.5*width*x for x in range(13)]
mouse_xs = [x+width for x in human_xs]
xtick_xs = [human_xs[i]+0.5*width for i in range(13)]

humanConserLncNums = []
mouseConcerLncNums = []
simianLncNums = []
primateLncNums = []
rodentLncNums = []
xlabels = []
for nums in numbers:
	humanConserLncNums.append(nums[3])
	mouseConcerLncNums.append(nums[5])
	simianLncNums.append(nums[1])
	primateLncNums.append(nums[2])
	rodentLncNums.append(nums[4])
	xlabels.append(nums[0].replace('_',' '))

plt.subplots(figsize=(5,3.5))
	
plt.bar(human_xs,
		primateLncNums,
		color='#F08080',#lightcoral
		width=width,
		label='Primate-specific'
		)
plt.bar(human_xs,
		simianLncNums,
		color='#8B0000',#darkred
		width=width,
		bottom=primateLncNums,
		label='Simian-specific'
		)
plt.bar(mouse_xs,
		rodentLncNums,
		color='#4169E1',#royalblue
		width=width,
		label='Rodent-specific'
		)
plt.yticks([0,300,600,900,1200],
		   [0,300,600,900,1200],fontsize='xx-small')
plt.xticks(xtick_xs, xlabels, fontsize='small', 
		   rotation=-55, verticalalignment='top', 
		   horizontalalignment='left')
plt.xlim(0.5-0.6*width, max(mouse_xs)+0.6*width)
plt.ylabel('Numbers of expressed lncRNAs', fontsize='small')
plt.legend(ncol=3, loc='lower left',fontsize='small',bbox_to_anchor=(-0.1, 1))
plt.tight_layout()
plt.savefig('fig/cancerExpressedLnc.pdf')
plt.close()