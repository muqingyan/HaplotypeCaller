import sys
#infile = open(sys.argv[1],'r')
#outfile = open(sys.argv[2],'w')
def vcfFilter(vcfin,vcfout):
    print("Now filtering VCF file %s...Variants with DP>=30 were recorded below.\n" % vcfin)
    infile = open(vcfin,'r+')
    outfile1 = open("%s.DP10.vcf" % vcfout,'w')
    outfile2 = open("%s.DP20.vcf" % vcfout,'w')
    outfile3 = open("%s.DP30.vcf" % vcfout,'w')
    for line in infile:
        if line.startswith("#"):
            outfile1.write(line)
            outfile2.write(line)
            outfile3.write(line)
        if not line.startswith("#"):
            line = line.strip()
            line = line.split('\t')
            alt = line[4]
            if not alt.__contains__(","):
                info = line[9]
                info = info.strip()
                info_list = info.split(':')
                # DP = info_list[2]
                AD = info_list[1]
                AD_list = AD.split(',')
                DP = info_list[2]
               # for item in AD_list:
               #     DP += int(item)
                AD_ref = AD_list[0]
                AD_alt = AD_list[1]
                try:
                    AF = float(AD_alt)/float(DP)
                except ZeroDivisionError:
                    pass
                else:
#                    if AF >= 0.2 and int(DP) >= 10 and int(AD_alt) >= 2:
                    if int(DP) >= 10:
                        outfile1.write("\t".join(line) + '\n')
                        if int(DP) >= 20:
                            outfile2.write("\t".join(line)+'\n')
                            if int(DP) >= 30:
                                outfile3.write("\t".join(line)+'\n')
                                print('\t'.join(line[0:5]) + '\t' + str(AD_alt)+';'+str(DP) + ';' + str(AF))
            else:
                print("A locus with multiple alternative alleles! Need to split first!\n")
   #             info = line[9]
   #             info = info.strip()
   #             info_list = info.split(':')
   #             # DP = info_list[2]
   #             AD = info_list[1]
   #             AD_list = AD.split(',')
   #             DP = 0
   #             for item in AD_list:
   #                 DP += int(item)
   #             AD_ref = AD_list[0]
   #             AD_list.pop(0)
   #             AD_alt = max(AD_list)
   #             # print(AD_alt,DP)
   #             try:
   #                 AF = float(AD_alt)/float(DP)
   #             except ZeroDivisionError:
   #                 pass
   #             else:
   #                 if AF >= 0.2 and int(DP) >= 10:
   #                     outfile.write("\t".join(line) + '\n')
   # 
