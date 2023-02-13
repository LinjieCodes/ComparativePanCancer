#!/usr/bin/env python

import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import decomposition as skldec
from sklearn import manifold

human_sample_num = {}
mouse_sample_num = {}

exp_dir = 'data/scaledGeneExp/'

#get human gene expression
print('Reading humam scaled gene expression')
human_exps = {}
files = [file for file in os.listdir(exp_dir) if 'human' in file]
for file in files:
	cancer = file.replace('human_', '')
	with open(exp_dir+file) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleNum = len(samples)
		human_sample_num[cancer] = sampleNum
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			for i in range(sampleNum):
				sample = samples[i]+'-'+cancer
				exp = expList[i]
				if sample not in human_exps:
					human_exps[sample] = {}
				human_exps[sample][gene] = exp
			
#get mouse gene expression
print('Reading mouse scaled gene expression')
mouse_exps = {}
files = [file for file in os.listdir(exp_dir) if 'mouse' in file]
for file in files:
	cancer = file.replace('mouse_', '')
	with open(exp_dir+file) as f:
		samples = f.readline().strip().split('\t')[1:]
		sampleNum = len(samples)
		mouse_sample_num[cancer] = sampleNum
		for line in f:
			cols = line.strip().split('\t')
			gene = cols[0]
			expList = [float(exp) for exp in cols[1:]]
			for i in range(sampleNum):
				sample = samples[i]+'-'+cancer
				exp = expList[i]
				if sample not in mouse_exps:
					mouse_exps[sample] = {}
				mouse_exps[sample][gene] = exp
	
#convert to dataframe
#sample row, gene column -- t-SNE required
print('df convert')
human_exps = pd.DataFrame(human_exps).T
mouse_exps = pd.DataFrame(mouse_exps).T

#colors for each cancer or tissue
colors = {'bladder_cancer':'#800080',#purple
		  'breast_cancer':'#BC8F8F',#rosybrown
		  'glioblastoma':'#8B0000',#drakred
		  'glioma':'#FAA460',#sandybrown
		  'kidney_cancer':'#8B4513',#saddlebrown
		  'leukemia':'#FFD700',#gold
		  'liver_cancer':'#FFFF00',#yellow
		  'lung_cancer':'#808000',#olive
		  'lymphoma':'#9ACD32',#yellowgreen
		  'ovarian_cancer':'#0000FF',#blue
		  'pancreatic_cancer':'#FF0000',#red
		  'prostate_cancer':'#00008B',#darkblue
		  'skin_melanoma':'#6A5ACD',#slateblue
		  'normal_bladder':'#7CFC00',#lawngreen
		  'normal_blood':'#008000',#green
		  'normal_brain':'#7FFFD4',#aquamarine
		  'normal_breast':'#2F4F4F',#darkslategrey
		  'normal_kidney':'#00FFFF',#cyan
		  'normal_liver':'#ADD8E6',#lightblue
		  'normal_lung':'#00BFFF',#deepskyblue
		  'normal_ovary':'#4682B4',#steelblue
		  'normal_pancreas':'#008080',#teal
		  'normal_prostate':'#EE82EE',#violet
		  'normal_skin':'#1E90FF'}#dodgerblue

#t-SNE for human
print('PCA for human')
pca_50 = skldec.PCA(n_components=50)
human_exps_pca50 = pca_50.fit_transform(human_exps)
print('t-SNE for human')
tsne = manifold.TSNE(n_jobs=10)
human_tsne = tsne.fit_transform(human_exps_pca50)
cancers = []
cancer_sample_index = {}
num = 0
sample_index = {}
for sample in human_exps.index:
	this_cancer = sample.split('-')[-1]
	if this_cancer not in cancer_sample_index:
		cancers.append(this_cancer)
		cancer_sample_index[this_cancer] = [num, num]
	else:
		cancer_sample_index[this_cancer][1] = num
	sample_index[num] = sample
	num += 1
plt.figure(figsize=(3, 3))
f=open('human_tSNE', 'w')
for cancer in cancers:
	start, end =cancer_sample_index[cancer]
	color = colors[cancer]
	plt.scatter(human_tsne[start:end+1, 0], human_tsne[start:end+1, 1],
				c=color, lw=0, alpha=0.5, label=cancer.replace('_', ' '),s=20)
	for i in range(start, end+1):
		f.write('\t'.join([str(human_tsne[i,0]), str(human_tsne[i,1]), sample_index[i]])+'\n')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.xticks([])
plt.yticks([])
plt.savefig('t-SNE-human.pdf')
plt.close()
f.close()	


#t-SNE for mouse
print('PCA for mouse')
pca_50 = skldec.PCA(n_components=50)
mouse_exps_pca50 = pca_50.fit_transform(mouse_exps)
print('t-SNE for mouse')
tsne = manifold.TSNE(n_jobs=10)
mouse_tsne = tsne.fit_transform(mouse_exps_pca50)
cancers = []
cancer_sample_index = {}
num = 0
sample_index = {}
for sample in mouse_exps.index:
	this_cancer = sample.split('-')[-1]
	if this_cancer not in cancer_sample_index:
		cancers.append(this_cancer)
		cancer_sample_index[this_cancer] = [num, num]
	else:
		cancer_sample_index[this_cancer][1] = num
	sample_index[num] = sample
	num += 1
plt.figure(figsize=(3, 3))
f=open('mouse_tSNE', 'w')
for cancer in cancers:
	start, end =cancer_sample_index[cancer]
	color = colors[cancer]
	plt.scatter(mouse_tsne[start:end+1, 0], mouse_tsne[start:end+1, 1],
				c=color, lw=0, alpha=0.5, label=cancer.replace('_', ' '))
	for i in range(start, end+1):
		f.write('\t'.join([str(mouse_tsne[i,0]), str(mouse_tsne[i,1]), sample_index[i]])+'\n')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.xticks([])
plt.yticks([])
plt.savefig('t-SNE-mouse.pdf')
plt.close()
f.close()


print('######')
print('Human sample numbers:')
for cancer in human_sample_num:
	print(cancer, human_sample_num[cancer])
print('######')
print('Mouse sample numbers:')
for cancer in mouse_sample_num:
	print(cancer, mouse_sample_num[cancer])