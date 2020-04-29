####
# to keep frameshift deletion, frameshift insertion,  nonsynonymous SNV
#stopgain, stoploss, splicing variants (LOF)
###

import sys
infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')
header = infile.readline()
for line in infile:
    line = line.strip()
    line = line.split('\t')
    func_refGene = line[6]
    exonicfunc_refgene = line[7]
    if func_refGene == "exonic":
        if exonicfunc_refgene.startswith("frameshift") or exonicfunc_refgene.__contains__("nonsynonymous") or exonicfunc_refgene.__contains__("stopgain") or exonicfunc_refgene.__contains__("stoploss"):
            outfile.write("\t".join(line)+"\n")

    elif func_refGene == "splicing":
        outfile.write("\t".join(line)+"\n")

    elif func_refGene.__contains__("exonic;splicing"):
        if exonicfunc_refgene.startswith("frameshift") or exonicfunc_refgene.__contains__("nonsynonymous") or exonicfunc_refgene.__contains__("stopgain") or exonicfunc_refgene.__contains__("stoploss"):
            outfile.write("\t".join(line)+"\n")

