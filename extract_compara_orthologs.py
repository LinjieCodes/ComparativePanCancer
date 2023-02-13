#!/usr/bin/env python

import multiprocessing
import subprocess
import re
import os

def extract_compara(ENSG):
	cmd_pattern = "wget -q --header='Content-type:text/x-orthoxml+xml' 'https://rest.ensembl.org/homology/id/%s?target_taxon=10090;type=orthologues'  -O compara_files/%s"		
	cmd = cmd_pattern % (ENSG, ENSG)
	ENSMUSG = None
	subprocess.Popen(cmd, shell=True).wait()
	if os.path.exists("compara_files/" + ENSG):
		with open("compara_files/" + ENSG) as f:
			file_content = f.read()
			if 'ortholog_one2one' in file_content:
				ENSMUSG_re = re.search('geneId="(ENSMUSG.+?)"', file_content)
				if ENSMUSG_re:
					ENSMUSG = ENSMUSG_re.group(1)
	return ENSG, ENSMUSG
	
human_genes = set()
with open("SRP007483") as f:
	f.readline()
	for line in f:
		human_genes.add(line.split('\t')[0])

if not os.path.exists("compara_files/"):
	os.mkdir("compara_files/")
	
pool = multiprocessing.Pool(processes = 10)
orthologs = []
for ENSG in human_genes:
	orthologs.append(pool.apply_async(extract_compara, (ENSG, )))
pool.close()
pool.join()
with open("human-mouse-orthologs-one2one", 'w') as f:
	for res in orthologs:
		ENSG, ENSMUSG = res.get()
		if ENSMUSG:
			f.write('\t'.join((ENSG, ENSMUSG)) + '\n')