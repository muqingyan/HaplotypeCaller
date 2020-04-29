#!/usr/bin/python

import sys

in_f = open(sys.argv[1], 'r')
out_f = open(sys.argv[2], 'w')
for line_1 in in_f:
    line_1_cl = line_1.strip().split('\t')
    eas = line_1_cl[30].split(',')[2]
    eas = eas[4::]
    if 'InterVar: Uncertain' in line_1_cl[15] or 'InterVar: Patho' in line_1_cl[15] or 'InterVar: Likely Patho' in line_1_cl[15]:
        if eas == '.' or float(eas) <= 0.005:
            out_f.write(line_1)
