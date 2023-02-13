#!/usr/bin/env python

import os
from urllib.request import urlretrieve
import subprocess

#directory to save the raw count data
count_dir = 'data/count/'

#read the study ID
studies = set()
with open('sample_list_outlier_unfiltered') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		studies.add(cols[0])
		
#read the url
urls = {}
with open('study-list-mus_musculus') as f:
	f.readline()
	for line in f:
		columns = line.split('\t')
		studyID = columns[0]
		count_url = columns[8]
		tpm_url = columns[7]
		urls[studyID] = (count_url, tpm_url)

#download the raw counts		
for studyID in studies:
	count_url, tpm_url = urls[studyID]
	if not os.path.exists(count_dir+studyID):
		try:
			urlretrieve(count_url, count_dir+studyID)
			print(studyID, " count download done!")
		except:
			print('ERROR:', count_url)