#!/usr/bin/python

filein = "../BJCH_gnomad2.0.af0.2dp10.debenign.intervar.withid.zygosity_variant_numofindiv_dedup"
fileout = open('../BJCH_debenign_variant_noSB.bed','w')
with open (filein,'r') as variant_f:
    for line in variant_f:
        line = line.strip()
        line_cl = line.split('\t')
        variant_cl = line_cl[0].split('_')
        chrome,start,end = variant_cl[0],variant_cl[1],variant_cl[2]
        start = str(int(start) -1)
        fileout.write(chrome + '\t' + start+'\t'+end+'\n')

