#!/bin/bash
# run this script in phaseres and supply fasta file of a population as a first argument
# output is in snap_results/population
# the main information I was interested in was "ps	pn	 ps/pn" this is header for pnps.txt
# the metric it actually calculates is synonymous and non-synonymous nucleotide diversity because it actually takes all pairs of sequences 
# present in a fasta file, calculate number of synonymous changes/ synonymous sites and non-synonymous changes/ non-synonymous sites 
# then make average per all comparison.
# the SNAP.pl is a bit broken as it does not output dn ds (those corrected rates) I was playing with the code a bit in order to get more decimals
# this script has been written 22.11.2018
# and it assumes that SNAP was downloaded and extracted to /home/vlkofly/programs/Snap_program/
# https://www.hiv.lanl.gov/repository/aids-db/PROGS/Snap/

# watch out for codon alignment
outdir=snap_results


if [ ! -d "${outdir}" ]; then mkdir $outdir;fi

f=$1
pop=${f/.fasta/}
#pad=1

if [ ! -d "${outdir}/${pop}" ]; then mkdir ${outdir}/${pop}; fi
#head $f
echo $pop

awk 'NR%2{printf "%s ",$0;next;}1' $f | tr -d ">" > ${pop}.snap.tmp
sed 's/\s[ATGC]\{1\}/ /g' ${pop}.snap.tmp > ${pop}.snap.inp

mv ${pop}.snap.inp ${outdir}/${pop}/

cd ${outdir}/${pop}

rm *summary*
rm *codon*
rm *bacground*

/home/vlkofly/programs/Snap_program/SNAP.pl ${pop}.snap.inp

grep -v "Average" summary* | awk '{sum_ps+=$9; sum_pn+=$10; T+=1}END{avg_ps=sum_ps/NR;avg_pn=sum_pn/NR;pnps=avg_pn/avg_ps; printf "%f\t%f\t%f\n", avg_ps, avg_pn, pnps}' > $pop.pnps.txt

