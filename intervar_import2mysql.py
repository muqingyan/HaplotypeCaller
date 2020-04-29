import pymysql
import sys

infile = sys.argv[1]
#outfile = open(sys.argv[2],'w')

db_table = "WES"
conn = pymysql.connect(host='172.16.156.98', port=3306, user='root', passwd='abc123++',db='WES')
cursor = conn.cursor()

with open(infile) as inf:
    for line in inf:
        if not line.startswith('#'):
            line_list = line.strip().split('\t')
            chr,start,end,ref,alt = line_list[0],int(line_list[1]),int(line_list[2]),line_list[3],line_list[4]
            cmd = "select * from intervar where chr='%s' and start=%s and end=%s and ref='%s' and alt='%s'" % (chr,start,end,ref,alt)
            try:
                cursor.execute(cmd)
                conn.commit()
                item = cursor.fetchone()
                if item is None:
                    cmd_2 = "insert into intervar values ('%s',%s,%s,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" % (line_list[0],int(line_list[1]),int(line_list[2]),line_list[3],line_list[4],line_list[5],line_list[6],line_list[7],line_list[8],line_list[9],line_list[10],line_list[11],line_list[12],line_list[13],line_list[14],line_list[15],line_list[16],line_list[17],line_list[18],line_list[19],line_list[20],line_list[21],line_list[22],line_list[23],line_list[24],line_list[25],line_list[26],line_list[27],line_list[28],line_list[29],line_list[30],line_list[31],line_list[32],line_list[33])
                    cursor.execute(cmd_2)
                    conn.commit()
            except Exception as err:
                conn.rollback()
                print(err)
conn.close()
