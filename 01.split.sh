#!/bin/bash
if [ $# != 1 ];then
	echo "Please use the program as: sh $0 <DIFF/INTS/NA/HDGC/>"
	exit 1;
fi
main="/nfs/home/Mandy/BJCH_190701_HGC"
master="${main}/result"
result="${main}/vsUNDEF/${1}"
for i in `cat ${main}/${1}sample.list`
	do cp ${master}/${i}*zygosity ${result}
done 
