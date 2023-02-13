#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from sklearn.linear_model import LinearRegression

#gene examples, Spearman rho and FDR (CIBERSORTx and quanTIseq)
genes = [('CXCL13', 'bladder_cancer', 'T cell CD8+', 0.65, 2.4e-65, 0.61, 9.9e-39),
		 ('CXCL10', 'lymphoma', 'Macrophage M1', 0.89, 3.7e-33, 0.73, 1.1e-16),
		 ('IL1RL1', 'prostate_cancer', 'Myeloid dendritic cell', 0.31, 9.9e-19, 0.34, 1.1e-22),
		 ('NOTCH1', 'prostate_cancer', 'Myeloid dendritic cell', 0.36, 4.4e-26, 0.39, 8.4e-31),
		 ('HLA-DRA', 'glioblastoma', 'T cell CD8+', 0.49, 1.7e-11, 0.24, 3.6e-3)]
				
human_lncs = [('ITGB2-AS1', 'liver_cancer', 'T cell CD8+', 0.65, 2.9e-43, 0.55, 8.9e-28),
			  ('MSC-AS1', 'liver_cancer', 'Macrophage M2', 0.58, 2.7e-32, 0.36, 2.3e-11),
			  ('UBL7-AS1', 'pancreatic_cancer', 'T cell CD8+', 0.32, 1.9e-4, 0.34, 3.7e-05)]
					 
mouse_lncs = [('D430020J02Rik', 'leukemia', 'Neutrophil', 0.46, 8.5e-06, 0.65, 4.1e-12)]

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

#plot for divergent genes
for gene, cancer, targetCell, rho1, fdr1, rho2, fdr2 in genes:
	#ENSG-symbol
	geneids = set()
	with open('divergence_result/'+cancer) as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			if cols[2] == gene:
				geneids.add(cols[0])
				geneids.add(cols[1])

	#read human gene exp
	human_exps = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		humanSamples = f.readline().strip().split('\t')[1:]
		samNum = len(humanSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneids:
				exps = [float(e) for e in cols[1:]]
				for i in range(samNum):
					human_exps[humanSamples[i]] = exps[i]

	#read mouse gene exp
	mouse_exps = {}
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		mouseSamples = f.readline().strip().split('\t')[1:]
		samNum = len(mouseSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneids:
				exps = [float(e) for e in cols[1:]]
				for i in range(samNum):
					mouse_exps[mouseSamples[i]] = exps[i]

	#read cell proportion
	CIBERSORT_result = {}
	quanTIseq_result = {}
	cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
			 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
			 'NK cell','Monocyte','Tregs']
	for cell in cells:
		CIBERSORT_result[cell] = {}
		quanTIseq_result[cell] = {}
	with open('immuneCell_comparison/TIMER2.0/infiltration_estimation_for_tcga.csv') as f:
		for line in f:
			cols = line.strip().split(',')
			thisSample = cols[0]
			if thisSample in humanSamples:
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

	humanX = []
	human_cbs = []
	human_qts = []
	mouseX = []
	mouse_cbs = []
	mouse_qts = []
	for s in human_exps:
		if s in CIBERSORT_result[targetCell]:
			humanX.append(human_exps[s])
			human_cbs.append(CIBERSORT_result[targetCell][s])
			human_qts.append(quanTIseq_result[targetCell][s])
	for s in mouse_exps:
		if s in CIBERSORT_result[targetCell]:
			mouseX.append(mouse_exps[s])
			mouse_cbs.append(CIBERSORT_result[targetCell][s])
			mouse_qts.append(quanTIseq_result[targetCell][s])
	#linear regression
	model_cbs = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_cbs.fit(np.array(humanX+mouseX).reshape(-1,1), human_cbs+mouse_cbs)
	
	model_qts = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_qts.fit(np.array(humanX+mouseX).reshape(-1,1), human_qts+mouse_qts)
	
	if targetCell == 'T cell CD8+':
		cell_abbr = 'CD8+ T cell'
	elif targetCell == 'T cell CD4+':
		cell_abbr = 'CD4+ T cell'
	elif targetCell == 'Myeloid dendritic cell':
		cell_abbr = 'Myeloid DC'
	elif targetCell == 'Macrophage M1':
		cell_abbr = 'M1 Macrophage'
	elif targetCell == 'Macrophage M2':
		cell_abbr = 'M2 Macrophage'
	else:
		cell_abbr = targetCell
	
	#plot scatter
	outFile = 'fig/%s_%s_CBS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(humanX, human_cbs, s=0.7, color='r')
	plt.scatter(mouseX, mouse_cbs, s=0.7, color='b')
	plt.plot([min(humanX+mouseX), max(humanX+mouseX)],
			 model_cbs.predict(np.array([min(humanX+mouseX),max(humanX+mouseX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), fontstyle='italic', fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho1), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr1), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')
	plt.ylabel(cell_abbr+' (CBS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()
	
	outFile = 'fig/%s_%s_QTS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(humanX, human_qts, s=0.7, color='r')
	plt.scatter(mouseX, mouse_qts, s=0.7, color='b')
	plt.plot([min(humanX+mouseX), max(humanX+mouseX)],
			 model_qts.predict(np.array([min(humanX+mouseX),max(humanX+mouseX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), fontstyle='italic', fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho2), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr2), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')	
	plt.ylabel(cell_abbr+' (QTS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()
	
###################
###################	
#plot for human lncs
for gene, cancer, targetCell, rho1, fdr1, rho2, fdr2 in human_lncs:
	#ENSG-symbol
	geneids = set()
	with open('lncrna_analysis/primateSpecific_lncRNA_summary') as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			if cols[1] == gene:
				geneids.add(cols[0])

	#read human gene exp
	human_exps = {}
	with open('data/scaledGeneExp/human_'+cancer) as f:
		humanSamples = f.readline().strip().split('\t')[1:]
		samNum = len(humanSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneids:
				exps = [float(e) for e in cols[1:]]
				for i in range(samNum):
					human_exps[humanSamples[i]] = exps[i]

	#read cell proportion
	CIBERSORT_result = {}
	quanTIseq_result = {}
	cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
			 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
			 'NK cell','Monocyte','Tregs']
	for cell in cells:
		CIBERSORT_result[cell] = {}
		quanTIseq_result[cell] = {}
	with open('immuneCell_comparison/TIMER2.0/infiltration_estimation_for_tcga.csv') as f:
		for line in f:
			cols = line.strip().split(',')
			thisSample = cols[0]
			if thisSample in humanSamples:
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
				
	humanX = []
	human_cbs = []
	human_qts = []
	for s in human_exps:
		if s in CIBERSORT_result[targetCell]:
			humanX.append(human_exps[s])
			human_cbs.append(CIBERSORT_result[targetCell][s])
			human_qts.append(quanTIseq_result[targetCell][s])
	#linear regression
	model_cbs = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_cbs.fit(np.array(humanX).reshape(-1,1), human_cbs)
	
	model_qts = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_qts.fit(np.array(humanX).reshape(-1,1), human_qts)
	
	if targetCell == 'T cell CD8+':
		cell_abbr = 'CD8+ T cell'
	elif targetCell == 'T cell CD4+':
		cell_abbr = 'CD4+ T cell'
	elif targetCell == 'Myeloid dendritic cell':
		cell_abbr = 'Myeloid DC'
	elif targetCell == 'Macrophage M1':
		cell_abbr = 'M1 Macrophage'
	elif targetCell == 'Macrophage M2':
		cell_abbr = 'M2 Macrophage'
	else:
		cell_abbr = targetCell
	
	#plot scatter
	outFile = 'fig/%s_%s_CBS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(humanX, human_cbs, s=0.7, color='r')
	plt.plot([min(humanX), max(humanX)],
			 model_cbs.predict(np.array([min(humanX),max(humanX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), 
	#		  fontstyle='italic', 
	#		  fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho1), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr1), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')
	plt.ylabel(cell_abbr+' (CBS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()
	
	outFile = 'fig/%s_%s_QTS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(humanX, human_qts, s=0.7, color='r')
	plt.plot([min(humanX), max(humanX)],
			 model_qts.predict(np.array([min(humanX),max(humanX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), 
	#		  fontstyle='italic', 
	#		  fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho2), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr2), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')	
	plt.ylabel(cell_abbr+' (QTS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()
	
	
######################
######################
#plot for mouse lncRNAs
for gene, cancer, targetCell, rho1, fdr1, rho2, fdr2 in mouse_lncs:
	#ENSG-symbol
	geneids = set()
	with open('lncrna_analysis/rodentSpecific_lncRNA_summary') as f:
		f.readline()
		for line in f:
			cols = line.split('\t')
			if cols[1] == gene:
				geneids.add(cols[0])

	#read mouse gene exp
	mouse_exps = {}
	with open('data/scaledGeneExp/mouse_'+cancer) as f:
		mouseSamples = f.readline().strip().split('\t')[1:]
		samNum = len(mouseSamples)
		for line in f:
			cols = line.strip().split('\t')
			if cols[0] in geneids:
				exps = [float(e) for e in cols[1:]]
				for i in range(samNum):
					mouse_exps[mouseSamples[i]] = exps[i]

	#read cell proportion
	CIBERSORT_result = {}
	quanTIseq_result = {}
	cells = ['B cell','T cell CD4+','T cell CD8+','Neutrophil',
			 'Macrophage M1','Macrophage M2','Myeloid dendritic cell',
			 'NK cell','Monocyte','Tregs']
	for cell in cells:
		CIBERSORT_result[cell] = {}
		quanTIseq_result[cell] = {}
				
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

	mouseX = []
	mouse_cbs = []
	mouse_qts = []
	for s in mouse_exps:
		if s in CIBERSORT_result[targetCell]:
			mouseX.append(mouse_exps[s])
			mouse_cbs.append(CIBERSORT_result[targetCell][s])
			mouse_qts.append(quanTIseq_result[targetCell][s])
	#linear regression
	model_cbs = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_cbs.fit(np.array(mouseX).reshape(-1,1), mouse_cbs)
	
	model_qts = LinearRegression(copy_X=True,
							 fit_intercept=True,
							 n_jobs=5,
							 normalize=False)
	model_qts.fit(np.array(mouseX).reshape(-1,1), mouse_qts)
	
	if targetCell == 'T cell CD8+':
		cell_abbr = 'CD8+ T cell'
	elif targetCell == 'T cell CD4+':
		cell_abbr = 'CD4+ T cell'
	elif targetCell == 'Myeloid dendritic cell':
		cell_abbr = 'Myeloid DC'
	elif targetCell == 'Macrophage M1':
		cell_abbr = 'M1 Macrophage'
	elif targetCell == 'Macrophage M2':
		cell_abbr = 'M2 Macrophage'
	else:
		cell_abbr = targetCell
	
	#plot scatter
	outFile = 'fig/%s_%s_CBS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(mouseX, mouse_cbs, s=0.7, color='b')
	plt.plot([min(mouseX), max(mouseX)],
			 model_cbs.predict(np.array([min(mouseX),max(mouseX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), fontstyle='italic', fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho1), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr1), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2-0.1
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')
	plt.ylabel(cell_abbr+' (CBS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()
	
	outFile = 'fig/%s_%s_QTS.pdf' % (gene, targetCell.replace(' ','_').replace('+',''))
	plt.figure(figsize=(2, 2))
	plt.scatter(mouseX, mouse_qts, s=0.7, color='b')
	plt.plot([min(mouseX), max(mouseX)],
			 model_qts.predict(np.array([min(mouseX),max(mouseX)]).reshape(-1,1)), 
			 color='k', linewidth=0.2)
	plt.xticks(fontsize=6)
	plt.yticks(fontsize=6)
	ax = plt.gca()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#plt.title('%s, %s' % (gene, cancer.replace('_', ' ')), fontstyle='italic', fontsize='small')
	title = '%s, %s' % (gene, cancer.replace('_', ' '))
	plt.text(0.1, 0.83, 'rho='+str(rho2), fontsize='small',transform=ax.transAxes)
	plt.text(0.1, 0.73, 'FDR='+str(fdr2), fontsize='small',transform=ax.transAxes)
	tx = 0-0.1*(len(title)-21)/2
	plt.text(tx, 1.1, title, fontstyle='italic', fontsize='small',transform=ax.transAxes)
	plt.xlabel('Gene NX',fontsize='small')	
	plt.ylabel(cell_abbr+' (QTS)',fontsize='small')
	plt.tight_layout()
	plt.savefig(outFile)
	plt.close()