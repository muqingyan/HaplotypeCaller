#!/usr/bin/python
# -*- coding: UTF-8 -*-
#Author: Muqing Yan; Date: 2019.04.19
# usage: this program is to check the zygosity of each variant reported as LOF or pathogenic by Intervar, files used are
#1. Intervar result (with sample ID which includes all the final variants, e.g.:115.DP10.processed_gnomAD2.1.0.001.patho.intervar.withid)
# 2. avinput file (e.g.: 115.DP10.processed_gnomAD.avinput,
# which contains the AC ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
# AC=1 -> het; AC=2 -> hom. Herein rendered by "het/hom" in the sixth column
# AC will generally relate to the zygosity, and it's counting the number of variant alleles. In most situations, if the call is heterozygous (0/1), then AC will be 1; whereas, if it's homozygous (1/1), AC will be 2.)


import sys,getopt

def main(argv):
    inter = ''
    avin = ''
    out = ''
    try:
        opts, args = getopt.getopt(argv,"hi:a:o:",["intervar=","avinput=","outfile="])
        if len(opts) < 1:
            print("Warning: to see the usage of the program: " + sys.argv[0] + ' -h')
            sys.exit(0)
    except getopt.GetoptError:
        print(sys.argv[0]+' -i <intervar file> -a <avinput file> -o <output file>')
        sys.exit(0)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0] + ' -i <intervar file> -a <avinput file> -o <output file>')
            sys.exit(0)
        elif opt in ('-i','--intervar'):
            inter = arg
        elif opt in ('-a', '--avinput'):
            avin = arg
        elif opt in ('-o', '--outfile'):
            out = arg

    intervar = open(inter,'r')
    avinput = open(avin,'r')
    outfile = open(out,'w')

    variants = []
    variants_lines = intervar.readlines()

    for line in variants_lines:
        line = line.strip()
        line_list = line.split('\t')
        sample_ID = line_list[0]
        variant = '_'.join(line_list[1:6])
        variants.append(variant)

    for line in avinput:
        line = line.strip()
        line_list = line.split('\t')
        variant_2 = '_'.join(line_list[0:5])
        if variant_2 in variants:
            zygosity = line_list[5]
            index = variants.index(variant_2)
            variant_line = variants_lines[index].split('\t')
            outfile.write('\t'.join(variant_line[0:7]) + '\t' + zygosity + '\t' + '\t'.join(variant_line[7::]))


if __name__ == '__main__':
    main(sys.argv[1:])

