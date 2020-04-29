#!/usr/bin/python
import sys
var_f = open(sys.argv[1],'r')
pile_f = open(sys.argv[2],'r')
sum_f = open(sys.argv[3],'w')

dp4carry = 10
dp4noncarry = 20

var_cl = var_f.readlines()
pile_cl = pile_f.readlines()

sum_f.write("Variant\t#ofvalidindivs\n")

for item in var_cl:
    numcarry = 0
    numnoncarry = 0
    valid_sum = 0
    variant = item.split('\t')
    variant_cl = variant[0].split('_')
    var_target = '_'.join(variant_cl[0:2])
    for line in pile_cl:
        line_cl = line.split('\t')
        var_indiv = '_'.join(line_cl[1:3])
        indiv,dp = line_cl[0],line_cl[4]
        if var_indiv == var_target:
            if indiv in item:
                #if int(dp) >= dp4carry:
                valid_sum += 1 
            elif int(dp) >= dp4noncarry:
                valid_sum += 1
    sum_f.write(variant[0]+'\t'+str(valid_sum)+'\n')

