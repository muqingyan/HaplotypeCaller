#!/usr/bin/python
import sys
main = "/nfs/home/Mandy/BJCH_redo_190510"
f_in = open("%s/BJCH_LOF.intervar.vargenepair" % main,'r')
g_f = open("%s/GC62gene.list" % main, 'r')
gene = []
for line in g_f:
    gene.append(line.strip()) 
#print(gene)
for line in f_in:
#    for item in gene
    if line.strip().split('\t')[1] in gene:
        #pass
        print(line.strip())
