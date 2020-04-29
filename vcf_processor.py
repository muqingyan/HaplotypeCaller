import sys
import os
import vcf_split
import vcf_AC_DP
def main():
    vcf_split.vcf_split(sys.argv[1],sys.argv[2])
    vcf_AC_DP.vcfFilter(sys.argv[2],sys.argv[3])
if __name__ == "__main__":
   main()
