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
		   
#read sample IDs
sample_cancer = {}		   
for cancer in cancers:
	with open('data/scaledGeneExp/human_'+cancer) as f:
		samples = f.readline().strip().split('\t')[1:]
		for sample in samples:
			sample_cancer[sample] = cancer
			
#read TIMER2.0 results			
results = {}
cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
		 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
		 'NK cell','Monocyte','Tregs']
for cancer in cancers:
	for cell in cells:
		results[cancer+'|'+cell] = {'CIBERSORT score':[[],[]],#human, mouse
									'quanTIseq score':[[],[]]}
with open('TIMER2.0/infiltration_estimation_for_tcga.csv') as f:
	for line in f:
		cols = line.strip().split(',')
		sample = cols[0]
		if sample in sample_cancer:
			cancer = sample_cancer[sample]
			values = [float(value) for value in cols[1:] if value !='NA']
			if cancer=='leukemia':#TIMER score is not available for leukemia
				values = ['NA','NA','NA','NA','NA','NA']+values
				
			results[cancer+'|B cell']['CIBERSORT score'][0].append(sum(values[28:31]))
			results[cancer+'|T cell CD4+']['CIBERSORT score'][0].append(sum(values[32:35]))
			results[cancer+'|T cell CD8+']['CIBERSORT score'][0].append(values[31])
			results[cancer+'|Neutrophil']['CIBERSORT score'][0].append(values[49])
			results[cancer+'|Macrophage M1']['CIBERSORT score'][0].append(values[42])
			results[cancer+'|Macrophage M2']['CIBERSORT score'][0].append(values[43])
			results[cancer+'|Myeloid dendritic cell']['CIBERSORT score'][0].append(sum(values[44:46]))
			results[cancer+'|Monocyte']['CIBERSORT score'][0].append(values[40])
			results[cancer+'|NK cell']['CIBERSORT score'][0].append(sum(values[38:40]))
			results[cancer+'|Tregs']['CIBERSORT score'][0].append(values[36])
			
			results[cancer+'|B cell']['quanTIseq score'][0].append(values[50])
			results[cancer+'|T cell CD4+']['quanTIseq score'][0].append(values[56])
			results[cancer+'|T cell CD8+']['quanTIseq score'][0].append(values[57])
			results[cancer+'|Neutrophil']['quanTIseq score'][0].append(values[54])
			results[cancer+'|Macrophage M1']['quanTIseq score'][0].append(values[51])
			results[cancer+'|Macrophage M2']['quanTIseq score'][0].append(values[52])
			results[cancer+'|Myeloid dendritic cell']['quanTIseq score'][0].append(values[59])
			results[cancer+'|Monocyte']['quanTIseq score'][0].append(values[53])
			results[cancer+'|NK cell']['quanTIseq score'][0].append(values[55])
			results[cancer+'|Tregs']['quanTIseq score'][0].append(values[58])
		
def sumLines(lines):
	sumList = [0]*len(lines[0].split(',')[1:])
	for line in lines:
		thisList = [float(value) for value in line.strip().split(',')[1:]]
		for index in range(len(thisList)):
			sumList[index]+= thisList[index]
	return sumList

for cancer in cancers:
	with open('TIMER2.0/mouse_%s.csv' % cancer) as f:
		lines = f.readlines()
		results[cancer+'|B cell']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0] in ('B cell naive_CIBERSORT-ABS','B cell memory_CIBERSORT-ABS','B cell plasma_CIBERSORT-ABS')])
		results[cancer+'|T cell CD4+']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0] in ('T cell CD4+ naive_CIBERSORT-ABS','T cell CD4+ memory resting_CIBERSORT-ABS','T cell CD4+ memory activated_CIBERSORT-ABS')])
		results[cancer+'|T cell CD8+']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='T cell CD8+_CIBERSORT-ABS'])
		results[cancer+'|Neutrophil']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Neutrophil_CIBERSORT-ABS'])
		results[cancer+'|Macrophage M1']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Macrophage M1_CIBERSORT-ABS'])
		results[cancer+'|Macrophage M2']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Macrophage M2_CIBERSORT-ABS'])
		results[cancer+'|Myeloid dendritic cell']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0] in ('Myeloid dendritic cell resting_CIBERSORT-ABS','Myeloid dendritic cell activated_CIBERSORT-ABS')])
		results[cancer+'|Monocyte']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Monocyte_CIBERSORT-ABS'])
		results[cancer+'|NK cell']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0] in ('NK cell resting_CIBERSORT-ABS','NK cell activated_CIBERSORT-ABS')])
		results[cancer+'|Tregs']['CIBERSORT score'][1] = sumLines([line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_CIBERSORT-ABS'])

		results[cancer+'|B cell']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='B cell_QUANTISEQ'])
		results[cancer+'|T cell CD4+']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='T cell CD4+ (non-regulatory)_QUANTISEQ'])
		results[cancer+'|T cell CD8+']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='T cell CD8+_QUANTISEQ'])
		results[cancer+'|Neutrophil']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Neutrophil_QUANTISEQ'])
		results[cancer+'|Macrophage M1']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Macrophage M1_QUANTISEQ'])
		results[cancer+'|Macrophage M2']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Macrophage M2_QUANTISEQ'])
		results[cancer+'|Myeloid dendritic cell']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Myeloid dendritic cell_QUANTISEQ'])
		results[cancer+'|Monocyte']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='Monocyte_QUANTISEQ'])
		results[cancer+'|NK cell']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='NK cell_QUANTISEQ'])
		results[cancer+'|Tregs']['quanTIseq score'][1] = sumLines([line for line in lines if line.split(',')[0]=='T cell regulatory (Tregs)_QUANTISEQ'])

#compare human and mouse

final_results = []
p_values = []
for cancer in cancers:
	for cell in cells:
		this_result = [cancer.replace('_',' '), cell]
		
		humanMean = np.mean(results[cancer+'|'+cell]['CIBERSORT score'][0])
		mouseMean = np.mean(results[cancer+'|'+cell]['CIBERSORT score'][1])
		logFC1 = np.log2((humanMean+0.0001)/(mouseMean+0.0001))	#add 0.0001 to avoid division by 0
		p1 = stats.mannwhitneyu(results[cancer+'|'+cell]['CIBERSORT score'][0], results[cancer+'|'+cell]['CIBERSORT score'][1])[1]
		this_result += [str(logFC1),str(p1),'NA',str(humanMean),str(mouseMean)]
		
		humanMean = np.mean(results[cancer+'|'+cell]['quanTIseq score'][0])
		mouseMean = np.mean(results[cancer+'|'+cell]['quanTIseq score'][1])
		logFC2 = np.log2((humanMean+0.0001)/(mouseMean+0.0001))	#add 0.0001 to avoid division by 0
		p2 = stats.mannwhitneyu(results[cancer+'|'+cell]['quanTIseq score'][0], results[cancer+'|'+cell]['quanTIseq score'][1])[1]
		this_result += [str(logFC2),str(p2),'NA',str(humanMean),str(mouseMean)]
		
		final_results.append(this_result)
		p_values.append(p1)
		p_values.append(p2)
		
#Benjamini-Hochberg p-value correction
p_values = np.asfarray(p_values)
by_descend = p_values.argsort()[::-1]
by_orig = by_descend.argsort()
steps = float(len(p_values)) / np.arange(len(p_values), 0, -1)
fdrs = np.minimum(1, np.minimum.accumulate(steps * p_values[by_descend]))
for i in range(len(final_results)):
	final_results[i][4] = str(fdrs[by_orig][i*2])
	final_results[i][9] = str(fdrs[by_orig][i*2+1])
	
#write result	
with open('immuneInfiltration_compare_result','w') as f:
	for this_result in final_results:
		if float(this_result[2]) * float(this_result[7])>0 and (float(this_result[4])<0.01 or float(this_result[9])<0.01):
			f.write('\t'.join(this_result)+'\n')