#!/usr/bin/env python

from scipy.spatial import distance
import numpy as np
import pandas as pd
import random
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

#calculate distances between samples
def cal_distance(df):
	dist_metric = distance.pdist(df, 'euclidean')
	dist_metric = distance.squareform(dist_metric)
	return dist_metric

#analysis of similarities	
def anosim(df, dist_metric, label1, label2, cancer, permu_time, outFile):
	index_list1 = []
	index_list2 = []
	index = 0
	for sample in df.index:
		if label1 in sample:
			index_list1.append(index)
		elif label2 in sample:
			index_list2.append(index)
		index += 1
	with open(outFile, 'a') as f:
		R, r_b, r_w, sample_num = cal_R(index_list1, index_list2, 
										dist_metric)
		p_val = permutation(index_list1, index_list2, dist_metric, 
							R, permu_time)
		f.write('\t'.join([cancer, str(R), 
						   '<'+str(p_val), str(int(round(r_b))), 
						   str(int(round(r_w))), str(sample_num), 
						   str(permu_time)])+'\n')

#calculate R values						   
def cal_R(index_list1, index_list2, dist_metric):
	dist_rank = []
	merged_index = []
	r_within = []
	r_between = []
	for i in index_list1:
		merged_index.append((i, 'g1'))
	for i in index_list2:
		merged_index.append((i, 'g2'))
	merged_index.sort()
	for i, label1 in merged_index:
		for j, label2 in merged_index:
			if i < j:
				if label1 == label2:
					group_flag = 'w' #within the group
				else:
					group_flag = 'b' #between the groups
				dist_rank.append((dist_metric[i][j], group_flag))
	dist_rank.sort()
	for rank in range(len(dist_rank)):
		dist, group_flag = dist_rank[rank]
		if group_flag == 'w':
			r_within.append(rank+1)
		else:
			r_between.append(rank+1)
	sample_num = len(merged_index)
	r_b_mean = np.mean(r_between)
	r_w_mean = np.mean(r_within)
	R = (r_b_mean - r_w_mean)/(0.25 * sample_num * (sample_num-1))
	return R, r_b_mean, r_w_mean, sample_num
	
#permutation
def permutation(index_list1, index_list2, dist_metric, R, permu_time):
	count = 0
	for t in range(permu_time):
		all_index = index_list1 + index_list2
		random.shuffle(all_index)
		new_index1 = all_index[:len(index_list1)]
		new_index2 = all_index[len(index_list1):]
		r = cal_R(new_index1, new_index2, dist_metric)[0]
		if r >= R:
			count += 1
	p_val = (count+1) / permu_time
	return p_val

#main fucntion	
def main():
	outFile = 'anosim_result'
	if os.path.exists(outFile):
		os.remove(outFile)
		
	permu_time = 10000
		
	#directory of NX data	
	nx_dir = 'data/scaledGeneExp/'
	
	#read human-mouse 1:1 orthologs
	orthologs = {}
	with open('gene_annotation/one-to-one-orthologs') as f:
		for line in f:
			ensg, ensm = line.strip().split('\t')
			orthologs[ensg] = ensg+'-'+ensm
			orthologs[ensm] = ensg+'-'+ensm
	
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
			   
	for cancer in cancers:
		human_exps = {}
		mouse_exps = {}
		with open(nx_dir+'human_'+cancer) as f:
			samples = [s+'-human' for s in f.readline().strip().split('\t')[1:]]
			sampleNum = len(samples)
			for line in f:
				cols = line.strip().split('\t')
				ensg = cols[0]
				if ensg in orthologs:
					ortho = orthologs[ensg]
					human_exps[ortho] = {}
					for i in range(sampleNum):
						human_exps[ortho][samples[i]] = float(cols[i+1])
		with open(nx_dir+'mouse_'+cancer) as f:
			samples = [s+'-mouse' for s in f.readline().strip().split('\t')[1:]]
			sampleNum = len(samples)
			for line in f:
				cols = line.strip().split('\t')
				ensm = cols[0]
				if ensm in orthologs:
					ortho = orthologs[ensm]
					mouse_exps[ortho] = {}
					for i in range(sampleNum):
						mouse_exps[ortho][samples[i]] = float(cols[i+1])
		
		#sample row, gene column 
		human_exps = pd.DataFrame(human_exps)
		mouse_exps = pd.DataFrame(mouse_exps)
		
		#merge two dataframes
		merged_df = pd.concat([human_exps, mouse_exps], axis=0).dropna(axis=1)
		
		#compute distance
		dist_metric = cal_distance(merged_df)

		#ANOSIM
		anosim(merged_df, dist_metric, 'human', 'mouse', cancer, permu_time, outFile)
		
if __name__ == '__main__':
	main()