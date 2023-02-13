#!/usr/bin/env python

import matplotlib.pyplot as plt

cladeLnc_ratios = []#######
simianLnc_ratios = []
primateLnc_ratios = []
rodentLnc_rations = []
n1 = 0
n2 = 0
n3 = 0
n4 = 0
with open('co_expression_network/coExpModule') as f:
	f.readline()
	for line in f:
		cols = line.split('\t')
		diverNum = int(cols[4])
		simianNum = int(cols[6])
		primateNum = int(cols[8])
		rodentNum = int(cols[10])
		
		cladeLnc_ratios.append((simianNum+primateNum+rodentNum)/diverNum)########
		simianLnc_ratios.append(simianNum/diverNum)
		primateLnc_ratios.append(primateNum/diverNum)
		rodentLnc_rations.append(rodentNum/diverNum)
		
		if (simianNum+primateNum+rodentNum)/diverNum>0.5:
			n1 += 1
		if (simianNum/diverNum>0.5 or primateNum/diverNum>0.5) and rodentNum/diverNum<0.1:
			n2 += 1
		if simianNum/diverNum<0.1 and primateNum/diverNum<0.1 and rodentNum/diverNum>0.5:
			n3 += 1
		if (simianNum+primateNum+rodentNum)/diverNum>0.8:
			n4 += 1
print(n1/len(simianLnc_ratios))
print(n2/len(simianLnc_ratios))
print(n3/len(simianLnc_ratios))
print(n4/len(simianLnc_ratios))

#histogram plotting
plt.figure(figsize=(2, 2.2))

width=0.6
xs = [1,2,3,4]

percentages = []

n2 = 0############
for r in cladeLnc_ratios:############
	if r>0.8:
		n2+=1
percentages.append(n2/len(cladeLnc_ratios))############
		
n2 = 0
for r in primateLnc_ratios:
	if r>0.8:
		n2+=1
percentages.append(n2/len(primateLnc_ratios))

n2 = 0
for r in simianLnc_ratios:
	if r>0.8:
		n2+=1
percentages.append(n2/len(simianLnc_ratios))

n2 = 0
for r in rodentLnc_rations:
	if r>0.8:
		n2+=1
percentages.append(n2/len(rodentLnc_rations))

xlabels = ['All clade-',############
		   'Primate-',
		   'Simian-',
		   'Rodent-']
plt.bar(xs,
		percentages,
		color='#1E90FF',#dodgerblue
		width=width,
		label='Ratio>0.8'
		)
		
for i in range(4):############
	p2 = percentages[i]
	plt.text(xs[i]-0.4*width, p2+0.01, str(int(p2*100))+'%',fontsize='x-small')
		
plt.yticks([0,0.25,0.5,0.75],
		   [0,25,50,75],fontsize='xx-small')
plt.xticks(xs, 
		   xlabels, 
		   fontsize='x-small',
		   rotation=-30,
		   #verticalalignment='top', 
		   #horizontalalignment='left'
		   )
#plt.xlim(0.5-0.6*width, max(xs2)+0.6*width)
plt.ylim(0,0.82)
plt.ylabel('Module percentage (%)', fontsize='small')
plt.xlabel('Clade-specific lncRNAs', fontsize='small')
plt.tight_layout()
plt.savefig('fig/percentages_module_enrichedWithLnc.pdf')
plt.close()