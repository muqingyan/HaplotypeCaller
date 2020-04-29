import sys
import re

def vcf_split(vcf_in,vcf_out):	
    print("Now processing VCF file %s...\n" % vcf_in)
    filein = open(vcf_in,'r')
    fileout = open(vcf_out,'w+')
    for line in filein:
        if line.startswith('#'):
            fileout.write(line)
        else:
            line = line.strip()
            line_list = line.split('\t')
            chr,pos,ref,alt = line_list[0],line_list[1],line_list[3],line_list[4]
            if not alt.__contains__(','):
                fileout.write(line+'\n')
            else:
                length = len(alt.split(','))
                alt_list = alt.split(',')
                AC_list = re.split(r'[=,]',line_list[7].split(';')[0])
                AF_list = re.split(r'[=,]',line_list[7].split(';')[1])
                GT,AD,DP,GQ,PL = line_list[9].split(':')[0], line_list[9].split(':')[1],line_list[9].split(':')[2],line_list[9].split(':')[3],line_list[9].split(':')[4]
                for i in range(0,length):
                    num = i*3
                    new_ac = 'AC='+ AC_list[i+1]
                    new_af = 'AF=' + AF_list[i+1]
                # new_GT = GT
                    new_AD = AD.split(',')[0] +','+AD.split(',')[i+1]
                    new_PL = ','.join(PL.split(',')[num:num+3])
                    new_AD_list = new_AD.split(',')
                    new_AD_list = map(int,new_AD_list)
                    new_DP = sum(new_AD_list)
                    new_line = '\t'.join(line_list[0:4])
                    new_line = new_line + '\t' + alt_list[i] + '\t'+'\t'.join(line_list[5:7])
                    new_line = new_line + '\t' + new_ac+';'+new_af+';'+';'.join(line_list[7].split(';')[2:]) + '\t'+line_list[8]
                    new_line = new_line + '\t' + GT+':'+new_AD+':'+str(new_DP)+':'+GQ+':'+new_PL
                    fileout.write(new_line+'\n')





