#!/usr/bin/python
#-*- coding: UTF-8 -*-
#Author: Muqing Yan; Date: 2019.04.25
import sys,getopt

def main(argv):
    variant = ''
    gene = ''
    out_f = ''
    try:
        opts,args = getopt.getopt(argv,"hv:g:o:",["variant=","gene=","outfile="])
        if len(opts) < 1:
            print("Warning: to see the usage of the program: " + sys.argv[0] + ' -h')
            sys.exit(0)
    except getopt.GetoptError:
        print('python ' + sys.argv[0] + ' -v <file_with_variants> -g <file_with_gene> -o <outfile>')
        sys.exit(0)
    for opt,arg in opts:
        if opt == '-h':
            print('python ' + sys.argv[0] + ' -v <file_with_variants> -g <file_with_gene> -o <outfile>')
            sys.exit(0)
        if opt in ('-v','--variant'):
            variant = arg
        if opt in ('-g','--gene'):
            gene = arg
        if opt in ('-o','--outfile'):
            out_f = arg
    
    variant_f = open(variant, 'r')
    gene_f = open(gene, 'r')
    out_file = open(out_f, 'w')

    genes = gene_f.readlines()
    for line in variant_f:
        var = line.strip().split()[0]
        for line2 in genes:
            if var in line2:
                out_file.write(var + '\t' + line2.split()[0] +'\t' +'\t'.join(line.strip().split()[1::])+'\n')

if __name__ == '__main__':
    main(sys.argv[1:])

