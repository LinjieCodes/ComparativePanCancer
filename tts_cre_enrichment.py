#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import random
import numpy as np
from scipy.stats import fisher_exact

#human and mouse TTS
tts_cre_human = 0
tts_noncre_human = 0
tts_cre_mouse = 0
tts_noncre_mouse = 0
tts_lengths = []
with open('TTS') as f:
	f.readline()
	for line in f:
		cols = line.strip('\r\n').split('\t')
		ttsList = cols[4].split(',')
		ttsNum = len(ttsList)
		for tts in ttsList:
			ttsStart = int(tts[tts.find(':')+1:tts.find('-')])
			ttsEnd = int(tts[tts.find('-')+1:])
			ttsLen = ttsEnd-ttsStart+1
			tts_lengths.append(ttsLen)
		if cols[0] == 'Human':
			if cols[5]:
				creNum = len(cols[5].split(','))				
				tts_cre_human += creNum
				tts_noncre_human += (ttsNum-creNum)
			else:
				tts_noncre_human += ttsNum
		else:
			if cols[5]:
				creNum = len(cols[5].split(','))				
				tts_cre_mouse += creNum
				tts_noncre_mouse += (ttsNum-creNum)
			else:
				tts_noncre_mouse += ttsNum
				
tts_mean_length = int(np.mean(tts_lengths))
print('TTS mean length:',tts_mean_length)
				
#human random region
cres = {}
chr_range = {}
with open('gene_annotation/GRCh38-ccREs.bed') as f:
	for line in f:
		chr, start, end = line.split('\t')[:3]
		start = int(start)
		end = int(end)
		if chr not in cres:
			cres[chr] = []
		cres[chr].append((start, end))
for chr in cres:
	cres[chr].sort()
	chr_range[chr] = (cres[chr][0][0],cres[chr][-1][1])
chrs = list(chr_range.keys())
random_cre_num = 0
random_nonCre_num=0
for i in range(10000):
	chr = random.choice(chrs)
	rang_mix,range_max = chr_range[chr]
	randomStart = random.randint(rang_mix,range_max)
	randomEnd = randomStart+tts_mean_length-1
	overlap = 0
	for creStart, creEnd in cres[chr]:
		if creStart>randomEnd:
			break
		elif creEnd<randomStart:
			continue
		else:
			overlapLen = min(randomEnd,creEnd) - max(randomStart,creStart) + 1
			if overlapLen>20:
				overlap=1
				break
			else:
				continue
	if overlap:
		random_cre_num +=1
	else:
		random_nonCre_num+=1
		
tts_cre_percent_human = tts_cre_human/(tts_cre_human+tts_noncre_human)
print(tts_cre_percent_human)
randome_cre_percent_human = random_cre_num/(random_cre_num+random_nonCre_num)
oddsr, p_human = fisher_exact([[tts_cre_human,random_cre_num],[tts_noncre_human,random_nonCre_num]], alternative='greater')

#plot human TTSs' CRE enrichment
plt.figure(figsize=(1.4, 2))
labels = ['DBSs', 'Random\nregions']
bar_width = 0.5
plt.bar([1, 1.8],
		[tts_cre_percent_human,randome_cre_percent_human],
		color='#FFA500',#orange
		width=bar_width, 
		label='CRE overlapped'
		)
plt.bar([1, 1.8],
		[1-tts_cre_percent_human,1-randome_cre_percent_human],
		color='#6495ED',#cornflowerblue
		width=bar_width,
		bottom=[tts_cre_percent_human,randome_cre_percent_human],
		label='No CRE'
		)
anno_human = 'P = '+str(p_human)
plt.xticks([1, 1.8],labels,rotation=-30,fontsize='x-small')
plt.yticks([0,0.5,1],[0,50,100],fontsize=4)
plt.text(1, 1.05, anno_human, fontsize='x-small', fontstyle='oblique')
plt.hlines(y=1.02, xmin=0.9, xmax=1.9, linewidth=0.7)
plt.ylabel('Percentage (%)', fontsize='x-small')
plt.xlim(0.6,2.2)
#plt.legend(ncol=1, loc='lower left',fontsize='xx-small',bbox_to_anchor=(-0.6, 1))
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('fig/TTS_CREenrich_human.pdf')
plt.close()


#mouse random region
cres = {}
chr_range = {}
with open('gene_annotation/mm10-ccREs.bed') as f:
	for line in f:
		chr, start, end = line.split('\t')[:3]
		start = int(start)
		end = int(end)
		if chr not in cres:
			cres[chr] = []
		cres[chr].append((start, end))
for chr in cres:
	cres[chr].sort()
	chr_range[chr] = (cres[chr][0][0],cres[chr][-1][1])
chrs = list(chr_range.keys())
random_cre_num = 0
random_nonCre_num=0
for i in range(10000):
	chr = random.choice(chrs)
	rang_mix,range_max = chr_range[chr]
	randomStart = random.randint(rang_mix,range_max)
	randomEnd = randomStart+tts_mean_length-1
	overlap = 0
	for creStart, creEnd in cres[chr]:
		if creStart>randomEnd:
			break
		elif creEnd<randomStart:
			continue
		else:
			overlapLen = min(randomEnd,creEnd) - max(randomStart,creStart) + 1
			if overlapLen>20:
				overlap=1
				break
			else:
				continue
	if overlap:
		random_cre_num +=1
	else:
		random_nonCre_num+=1
		
tts_cre_percent_mouse = tts_cre_mouse/(tts_cre_mouse+tts_noncre_mouse)
print(tts_cre_percent_mouse)
randome_cre_percent_mouse = random_cre_num/(random_cre_num+random_nonCre_num)
oddsr, p_mouse = fisher_exact([[tts_cre_mouse,random_cre_num],[tts_noncre_mouse,random_nonCre_num]], alternative='greater')

#plot mouse TTSs' CRE enrichment
plt.figure(figsize=(1.4, 2))
labels = ['DBSs', 'Random\nregions']
bar_width = 0.5
plt.bar([1, 1.8],
		[tts_cre_percent_mouse,randome_cre_percent_mouse],
		color='#FFA500',#orange
		width=bar_width, 
		label='CRE overlapped'
		)
plt.bar([1, 1.8],
		[1-tts_cre_percent_mouse,1-randome_cre_percent_mouse],
		color='#6495ED',#cornflowerblue
		width=bar_width,
		bottom=[tts_cre_percent_mouse,randome_cre_percent_mouse],
		label='No CRE'
		)
anno_mouse = 'P = '+str(p_mouse)
plt.xticks([1, 1.8],labels,rotation=-30,fontsize='x-small')
plt.yticks([0,0.5,1],[0,50,100],fontsize=4)
plt.text(1, 1.05, anno_mouse, fontsize='x-small', fontstyle='oblique')
plt.hlines(y=1.02, xmin=0.9, xmax=1.9, linewidth=0.7)
plt.ylabel('Percentage (%)', fontsize='x-small')
plt.xlim(0.6,2.2)
#plt.legend(ncol=1, loc='lower left',fontsize='xx-small',bbox_to_anchor=(-0.6, 1))
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('fig/TTS_CREenrich_mouse.pdf')
plt.close()