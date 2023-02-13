#!/usr/bin/env python

import matplotlib.pyplot as plt

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
		   
#gene number counting
geneNums = {}
for cancer in cancers:
	geneNums[cancer] = {'divergent':0,
						'conserved':0,
						'other':0,
						'total':0}
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.strip().split('\t')
			if abs(float(cols[4]))>1 and float(cols[-1])<0.05:
				geneNums[cancer]['divergent'] += 1
			elif abs(float(cols[4]))<0.5:
				geneNums[cancer]['conserved'] += 1
			else:
				geneNums[cancer]['other'] += 1
			geneNums[cancer]['total'] += 1
				
#plotting
plt.subplots(figsize=(6,3.5))

width=0.85
x = list(range(13))
xlabels = []

divergent_percents = []
conserved_percents = []
other_percents = []
for cancer in cancers:
	xlabels.append(cancer.replace('_',' '))
	divergent_percents.append(geneNums[cancer]['divergent']/geneNums[cancer]['total'])
	conserved_percents.append(geneNums[cancer]['conserved']/geneNums[cancer]['total'])
	other_percents.append(geneNums[cancer]['other']/geneNums[cancer]['total'])
	
plt.bar(x,
		conserved_percents,
		color='#6A5ACD',#slateblue
		width=width,
		label='Conserved'
		)
plt.bar(x,
		other_percents,
		color='#FFA500',#orange
		width=width,
		bottom=conserved_percents,
		label='Intermediate'
		)
plt.bar(x,
		divergent_percents,
		color='#CD5C5C',#indianred
		width=width,
		bottom=[conserved_percents[i]+other_percents[i] for i in range(13)],
		label='Divergent'
		)
		
for i in range(13):
	cancer = cancers[i]
	y = 0.4*conserved_percents[i]
	num = geneNums[cancer]['conserved']
	if num>1000:
		plt.text(i-0.35*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	elif num<100:
		plt.text(i-0.2*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	else:
		plt.text(i-0.28*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	
	y = conserved_percents[i]+0.4*other_percents[i]
	num = geneNums[cancer]['other']
	if num>1000:
		plt.text(i-0.35*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	elif num<100:
		plt.text(i-0.2*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	else:
		plt.text(i-0.28*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	
	y = conserved_percents[i]+other_percents[i]+0.37*divergent_percents[i]
	num = geneNums[cancer]['divergent']
	if num>1000:
		plt.text(i-0.35*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	elif num<100:
		plt.text(i-0.2*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
	else:
		plt.text(i-0.28*width, y, str(num), fontsize='xx-small', fontstyle='oblique')
		
plt.yticks([0,0.25,0.5,0.75,1],
		   [0,25,50,75,100],fontsize='xx-small')
plt.xticks(x, xlabels, fontsize='small', 
		   rotation=-55, verticalalignment='top', 
		   horizontalalignment='left')
plt.ylabel('Gene percentage (%)', fontsize='small')
plt.legend(ncol=3, loc='center',fontsize='small',bbox_to_anchor=(0.5, 1.08))
plt.xlim(-0.6,12.6)
plt.tight_layout()
plt.savefig('fig/divergentGene_percentage.pdf')
plt.close()