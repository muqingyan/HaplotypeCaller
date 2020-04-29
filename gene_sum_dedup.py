#_*_coding:utf-8_*_
import sys

filein = open(sys.argv[1],'r')
variant_numofindiv = open("%s_variant_numofindiv_dedup" % sys.argv[1],'w')
indiv_numofgenes = open("%s_indiv_numofgenes_dedup" % sys.argv[1],'w')
gene_numofindivs  = open("%s_gene_numofindivs_dedup" % sys.argv[1],'w')
gene_numofvariants = open("%s_gene_numofvariants_dedup" % sys.argv[1],'w')
indiv_numofvariants = open("%s_indiv_numofvariants_dedup" % sys.argv[1],'w')
#genelist = open(sys.argv[7],'r')

#gene_list = []
#for line in genelist:
#    line = line.strip()
#    gene_list.append(line)

#header = filein.readline()
variant_dict = {}
gene_dict = {}
gene_variant_dict = {}
sample_variant_dict = {}
sample_gene_dict = {}
for line in filein:
    if not line.startswith('#'):
        line = line.strip()
        line_list = line.split('\t')
        key_variant = '_'.join(line_list[1:6])
        value_sample = line_list[0]
        key_gene = line_list[6]
        if key_variant not in variant_dict.keys():
            variant_dict[key_variant] = value_sample
        else:
            if value_sample not in variant_dict[key_variant]:
                variant_dict[key_variant] = variant_dict[key_variant] + ',' + value_sample
        if key_gene not in gene_dict.keys():
            gene_dict[key_gene] = value_sample
        else:
            if value_sample not in gene_dict[key_gene]:
                gene_dict[key_gene] = gene_dict[key_gene] + ',' + value_sample
        if key_gene not in gene_variant_dict.keys():
            gene_variant_dict[key_gene] = key_variant
        else:
            if key_variant not in gene_variant_dict[key_gene]:
                gene_variant_dict[key_gene] = gene_variant_dict[key_gene] + ',' + key_variant
        if value_sample not in sample_gene_dict:
            sample_gene_dict[value_sample]  = key_gene
        else:
            if key_gene not in sample_gene_dict[value_sample]:
                sample_gene_dict[value_sample] = sample_gene_dict[value_sample] + ',' + key_gene
        if value_sample not in sample_variant_dict:
            sample_variant_dict[value_sample] = key_variant
        else:
            if key_variant not in sample_variant_dict[value_sample]:
                sample_variant_dict[value_sample] = sample_variant_dict[value_sample] + ',' + key_variant
for key,value in variant_dict.items():
    length = len(value.split(','))
    variant_numofindiv.write(key + '\t' + value +'\t'+ str(length) +'\n')
for key,value in gene_dict.items():
    length = len(value.split(','))
    gene_numofindivs.write(key + '\t' + value +'\t'+ str(length)+'\n')
for key,value in gene_variant_dict.items():
    length = len(value.split(','))
    gene_numofvariants.write(key + '\t' + value +'\t'+ str(length)+'\n')
for key,value in sample_gene_dict.items():
    length = len(value.split(','))
    indiv_numofgenes.write(key + '\t' + value +'\t'+ str(length)+'\n')
for key,value in sample_variant_dict.items():
    length = len(value.split(','))
    indiv_numofvariants.write(key + '\t' + value +'\t'+ str(length)+'\n')
