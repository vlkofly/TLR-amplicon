#!/bin/bash

for f in *fasta
do
sed 's/>//g' $f | paste -d "\t" - - > ${f/fasta/sn}
done
