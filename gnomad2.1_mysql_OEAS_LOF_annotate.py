import pymysql
import sys

infile = sys.argv[1]
outfile = open(sys.argv[2],'w')

db_table = "gnomad2_1_otherEAST_LOF"
conn = pymysql.connect(host='172.16.156.98', port=3306, user='root', passwd='abc123++',db='GnomAD2_1')
cursor = conn.cursor()

with open(infile) as inf:
    for line in inf:
        line_list = line.strip().split('\t')
        chr,start,end,ref,alt = line_list[0],int(line_list[1]),int(line_list[2]),line_list[3],line_list[4]    
        cmd = "select otherinfo from %s where chr='%s' and start=%d and end=%d and ref='%s' and alt='%s'" % (db_table,chr,start,end,ref,alt)
        cursor.execute(cmd)
        AF_oea = cursor.fetchone()
        if AF_oea is None:
            af_oea = "." 
        else:
            af_oea = AF_oea[0].split('=')[1] 
        outfile.write('\t'.join(line_list[:5])+'\t'+af_oea+'\n') 
conn.close()
