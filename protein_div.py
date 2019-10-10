#!/usr/bin/python
"""This script takes a fasta AMK file created from geneious as first argument and parses it
The fasta translation is done on all DNA haplotypes, therefore non-unique amk sequences are present"""
# to do: make the translation automatic
import sys
import difflib


# variables

# parse user input
try:
    amkf = sys.argv[1]
except:
    print(__doc__)
    sys.exit(1)

# read input file
amf = open(amkf,'r')
aml = amf.readlines()
haploname = aml[0::2] # select haplotype names
amkseq = aml[1::2] # select amk sequences
uniqamk=set(amkseq) # by using set extract unique sequences
uniqamk = list(map(lambda s: s.strip("\n"), uniqamk)) # strip all the sequences of newline
#print(uniqamk)

# dictionary of protein haplotypes as key and haplonames as value
d = {}
#d.fromkeys(uniqamk)
for h,s in zip (haploname,amkseq):
    h = h.strip("\n").strip(">")
    #print(h, s)
    s = s.strip("\n")
    if s in d:
        d[s].append(h)
    else:
        d[s] = [h]

# export results
# segregating sites, get them from geneious after generating unique file, I improved the script it finds the positions
# by itself
# the position of segregating site starts from 0 from the beginning of my amplicon
# so later I should insert it into the whole protein and change the position
#div_pos = [20,21,23,39,43,46,50,58,73,164,172,222,225,242,245,270]

# loop through haplotypes and write out fasta files and tables:

with open('prot_table.tsv', 'w') as t: # table with  protein haplotype code, corresponding dna haplotypes and amk seq
    with open('prot_unique.fasta', 'w') as funiq: # fasta with the unique protein haplotypes
        with open('prot_segregating.fasta', 'w') as fdiv: # fasta with unique protein haplotypes limited to segregating sites, used as import of quantiprot
            with open('prot_segregating.csv', 'w') as p: # table with position of seg. site and amk, same format as the quantiprot metrics
                x = 0
                div_p = map(lambda z: [i for i in range(len(s)) if s[i] != z[i]], d)  # get the seg. sites in each comparison
                div_pos = list(set(z for y in div_p for z in y)) # unique seg. sites
                div_pos = sorted(div_pos) # sorted seg. sites
                #print(div_pu)
                t.write('prot_code\tdna_code\tprot_seq\n') # header of prot_table
                posstr = ','.join(map(str, div_pos))  # header of possitions
                p.write(',' + posstr + '\n')
                for seq in d: # loop through protein haplotypes
                    x += 1 # name of haplotype is just numeric
                    hap = d[seq] # dna haplotypes
                    hapstr = '-'.join(hap) # formating for export
                    hapstr_r = 'c('+','.join(hap)+')'
                    fname = '>'+str(x)+'_dnah'+hapstr+'\n'
                    funiq.write(fname)
                    funiq.write(seq+'\n')
                    t.write(str(x)+'\t'+hapstr_r+'\t'+seq+'\n')
                    fname_div = '>div' + str(x) + '\n'
                    ss = [seq[i] for i in div_pos]  # get AMK of segregating sites only
                    sj = ''.join(ss) # formating of AMK
                    sc = ','.join(ss)
                    fdiv.write(fname_div)
                    fdiv.write(sj + '\n')
                    p.write(str(x) + ',' + sc + '\n')





###quantiprot analysis, it will write out tables with AMK properties for segregating sites
# in order to see what are the physiochemical differences between the haplotypes
from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import SequenceSet
from quantiprot.utils.sequence import subset, columns
from quantiprot.utils.feature import Feature, FeatureSet

# Conversions-related imports:
from quantiprot.utils.mapping import simplify
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy, get_aa2volume,get_aa2mj
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.metrics.basic import identity
import numpy as np


fapr = load_fasta_file('prot_segregating.fasta') #load fasta
fs = FeatureSet("myTLRset")
fs.add(get_aa2charge())
fs.add(get_aa2volume())
fs.add(get_aa2mj())
fs.add(get_aa2hydropathy())

convfapr = fs(fapr)
metrics = ["formal_charge","volume","miyazawa-jernigan","hydropathy"] # which metrics of AMK to generate

#print convfapr
for m in metrics:
    outf = open(m+".tsv","w")
    with outf as f:
        h= np.matrix(columns(convfapr, feature=m, transpose=True))
        f.write(posstr+'\n')
        np.savetxt(f, h, delimiter=',')

