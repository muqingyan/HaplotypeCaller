#!/usr/bin/python
main = '/nfs/home/Mandy/BJCH_redo_190510'
var_f = open('%s/BJCH_LOF.intervar.varpairgeneindiv' % main, 'r')
intervar = open('%s/BJCH_LOF.intervar' % main, 'r')
#keys: 'gene', 'intervar','indiv'
var_dict = {'Variant(chromosome:start-end:ref>alt)':{'gene':'Gene','intervar':'Intervar','indiv':'Number of carriers'}}
intervar_cl = intervar.readlines()
for line in var_f:
    line_cl = line.strip().split('\t')
    variant,gene,numofindiv = line_cl[0],line_cl[1],line_cl[3]
    variant_cl = variant.split('_')
    variant = variant_cl[0]+':'+variant_cl[1]+'-'+variant_cl[2]+':'+variant_cl[3]+'>'+variant_cl[4]
    var_dict[variant] = {}
    var_dict[variant]['gene'] = gene
    var_dict[variant]['indiv'] = numofindiv
#    print(variant)
    for line_2 in intervar_cl:
        line_2cl = line_2.strip().split('\t')
        variant_2 = line_2cl[1]+':'+line_2cl[2]+'-'+line_2cl[3]+':'+line_2cl[4]+'>'+line_2cl[5]
        varinter = line_2cl[14]
        if variant_2 == variant:
#            print('===='+variant_2+'===')
            var_dict[variant]['intervar'] = varinter
#            print(var_dict[variant]['intervar'])
            break
#print(var_dict)
for key in var_dict.keys():
    print(var_dict[key]['gene']+'\t'+key+'\t'+var_dict[key]['indiv']+'\t'+var_dict[key]['intervar'])
        
