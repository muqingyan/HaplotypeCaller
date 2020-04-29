import pymysql
import sys

infile = sys.argv[1]
#outfile = open(sys.argv[2],'w')

db_table = "BJCH"
conn = pymysql.connect(host='172.16.156.98', port=3306, user='root', passwd='abc123++',db='BJCH')
cursor = conn.cursor()

with open(infile) as inf:
    for line in inf:
        if not line.startswith('#'):
            line_list = line.strip().split('\t')
            sample_id,chr,start,end,ref,alt,rs,qual,filter,info,format = line_list[0],line_list[1],int(line_list[2]),int(line_list[3]),line_list[4],line_list[5],line_list[6],line_list[7],line_list[8],line_list[9],line_list[11]
            cmd_1 = "select * from vcf where sample_id='%s' and chr='%s'and start=%s and end=%s and ref='%s' and alt='%s'" % (sample_id,chr,start,end,ref,alt)
            try:
                cursor.execute(cmd_1)
                conn.commit()
                item = cursor.fetchone()
                if item is None:
                    cmd_2 = "insert into vcf values ('%s','%s',%s,%s,'%s','%s','%s','%s','%s','%s','%s')" % (sample_id,chr,start,end,ref,alt,rs,qual,filter,info,format)
                    cursor.execute(cmd_2)
                    conn.commit()
            except Exception as err:
                conn.rollback()
                print(err)
conn.close()
