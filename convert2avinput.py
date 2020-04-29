#-*-coding:utf-8-*-
#!/usr/bin/python
#Author: Muqing Yan; Date: 2019-07-03
#This program is to help convert VCF files to avinput files to assist the downstream procedures.
import sys
import getopt

def adjustStartEndRefAlt(newstart, newend, newref, newalt):
    while newref[-1] == newalt[-1]:
        newref = newref[0:(len(newref)-1)]
        newalt = newalt[0:(len(newalt)-1)]
        newend -= 1
        if newref == '':
            newref = '-'
            newstart -= 1
            break
        if newalt == '':
            newal ='-';
            break

    while newref[0:1] == newalt[0:1]:
        newref = newref[1::]
        newalt = newalt[1::]
        newstart += 1
        if newref == '':
            newref = '-'
            newstart -= 1
            break
        if newalt == '':
            newalt = ''
            break
    return newstart, newend, newref, newalt

def convert(vcf_in, av_out):
    print('Processing sample...')
    vcf = open(vcf_in, 'r')
    av = open(av_out, 'w')
    # av.write('test')
    for line_1 in vcf:
        if line_1.startswith('#CHROM'):
            sample = line_1.strip().split('\t')[-1]
            print(sample)
        elif not line_1.startswith('#'):
            line_1_cl = line_1.strip().split('\t')
            chr, start, rs, ref, alt = line_1_cl[0], int(line_1_cl[1]), line_1_cl[2], line_1_cl[3], line_1_cl[4]
            if len(ref) == 1 and len(alt) ==1:
                newstart, newend = start, start+len(ref)-1
                newref, newalt = ref, alt
            elif len(ref) > len(alt):
                head = ref[0:len(alt)]
                if head == alt:
                    newstart, newend = start + len(head), start + len(ref) -1
                    newref, newalt = ref[len(alt)::], '-'
                else:
                    newstart, newend = start, start+len(ref) -1
                    newref, newalt = ref, alt
                newstart, newend, newref, newalt = adjustStartEndRefAlt(newstart, newend, newref, newalt)
            elif len(ref) < len(alt):
                head = alt[0:len(ref)]
                if head == ref:
                    newstart, newend = start+len(ref)-1, start+len(ref)-1
                    newref, newalt = '-', alt[len(ref)::]
                else:
                    newstart, newend = start, start+len(ref) -1
                    newref, newalt = ref, alt
                newstart, newend, newref, newalt = adjustStartEndRefAlt(newstart, newend, newref, newalt)
            else:
                head = ref[0:(len(ref)-1)]
                if alt.startswith(head):
                    newstart, newend = start + len(ref) -1, start+len(ref) -1
                    chopped_ref = ref
                    newref, newalt = chopped_ref[0:(len(chopped_ref) -1)], alt[0:(len(alt) -1)]
                else:
                    newstart, newend = start, start + len(ref)-1
                    newref, newalt = ref, alt
                newstart, newend, newref, newalt = adjustStartEndRefAlt(newstart, newend, newref, newalt)
            av.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sample, chr, newstart, newend, newref, newalt, rs,'\t'.join(line_1_cl[5::])))
def main(argv):
    vcf_in = ''
    av_out = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["vcf=","avout="])
        if len(opts) < 1:
            print("Warning: to see the usage of the program: python " + sys.argv[0] + ' -h')
            sys.exit(0)
    except getopt.GetoptError:
        print("Please use the program as: python " + sys.argv[0] + ' -i <.vcf> -o <.avinput>')
        sys.exit(0)
    for opt, arg in opts:
        if opt == '-h':
            print("Please use the program as: python " + sys.argv[0] + ' -i <.vcf> -o <.avinput>')
            sys.exit(0)
        elif opt in ('-i', '--vcf'):
            vcf_in = arg
        elif opt in ('-o', '--avout'):
            av_out = arg

    convert(vcf_in, av_out)

if __name__ == '__main__':
    main(sys.argv[1:])
