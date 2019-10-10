#! /bin/bash
# COLLAPSING
for i in `cat samplepass.txt`
do
echo $i

cat $i'.fastq' | fastx_collapser -v -Q 33 | sed -e 's#>#>tmpcontig_#g' -e 's#-#_freq_#g'  > '../collapsed'$i'.fasta'
done | tee  ../collapsed/collapse1.log 

