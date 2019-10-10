#! /bin/bash
#this script takes pair end sequences and assembles them 
date +"%T"
for f in *paired/*P.gz
do
if [[ $f == *1P.gz ]] 
then
r=${f/1P/2P} 
o=${f/_1fastq_1P.gz/}
o=${o/paired\//}
#echo $o
/home/vlkofly/programs/pear-0.9.10-bin-64/pear-0.9.10-bin-64 \
-v 10 -n 300 -q 30 -f $f -r $r -o '../assembled/'$o #-v 100bp overlap -n 300bp minimal size of assembly
fi 
done | tee  ../assembled/assembly2.log 
