import sys
infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')
#header = infile.readline()
for line in infile:
    line = line.strip()
    line_list = line.split('\t')
    EAS_AF = line_list[5]
    if EAS_AF == ".":
        outfile.write("\t".join(line_list[0:5]) + '\t' + EAS_AF +'\n')
    elif round(float(EAS_AF),3) <= 0.001:
        outfile.write("\t".join(line_list[0:5]) + '\t' + str(EAS_AF) +'\n')
