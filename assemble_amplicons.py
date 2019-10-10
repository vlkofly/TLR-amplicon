#! /usr/bin/python3
import os
import itertools
import glob
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
#from Bio.Align.Applications import MuscleCommandline
#cline = MuscleCommandline(input='example_protein.fasta', 
#                                            out='example_protein_aligned.fasta', 
#                                             verbose=True)

import difflib

# run this script from folder true
# script takes the amplicons per individual and assembles them per TLR
# make list of samples:  ls trueMT* | cut -d. -f 1 | uniq | sort > samplelist.txt
samlist = open('samplelist.txt','r')

total_N_samam = 0
assembled_N_samam = 0

for nn, s in enumerate(samlist): # loop throgh files and the two TLR
 s = s.strip("\n")
 #s = s[0:1]
 for tlr in [".-TLR15",".-TLR4"]:
  filestr = s+tlr+'*.fasta'
  print ('processing: '+s+tlr)
  files = glob.glob(filestr)
  print ('found '+str(len (files)))#'files: '+str(files))
  total_N_samam += 1
  seqs = []
  Nheterozygots = 0 # set iterators for how many homo/heterozygots are there in the files
  Nhomozygots = 0 
  for f in files: # open each file with the amplicon sequences and write the sequences to variable

   with open(f,'r') as fi:
    lines = fi.readlines()
    linelen = len(lines)
    if linelen == 4:# append sequences if there are two of them in the file, meaning the amplicon is heterozygous
     print('genotype '+f+' heterozygot')
     Nheterozygots += 1
     al1 = (lines[1].strip('\n'))
     al2 = (lines[3])
     #print(al1[290:465]+'\n'+al2[290:465])
     print((al1)+'\n'+(al2))
     #seqs.append(lines[1].strip('\n'))
     #seqs.append(lines[3])

     hetdif = difflib.ndiff(al1,al2)     #check where are the differences in between the heterozygots
     difs = [(p,n) for n,p in enumerate(hetdif) if p[0] in ('+','-')]
     print (difs)
    else:
     print('genotype '+f+' homozygot')
     Nhomozygots += 1
     seqs.append(lines[1])
  print(str(Nheterozygots)+' heterozygous amplicons; '+str(Nhomozygots)+' homozygous amplicons')
  nseqs = (len(seqs))
  if (nseqs == 3 and Nheterozygots == 0): # do assembly for homozygots in all three amplicons
   print ("processing homozygous only Assembly:"+s+tlr)
   asem = Assembly((Dseqrecord(seqs[0]),Dseqrecord(seqs[1]),Dseqrecord(seqs[2])), only_terminal_overlaps = True, limit = 14)
   if asem.linear_products == []:
    print (s+tlr+' assembly failed\n\n') # if I do it for only homozygous sequences, than no fails
   else:

    assembled_N_samam += 1
    aswat = asem.linear_products[0].seq.watson
    aswatlen = len(aswat)
    print(s+tlr+'assembly successful length: '+str(aswatlen)+'\n\n')
    with open('../finalassem/asemseq'+s+tlr+'.fasta','w') as asfile:
     asfile.write(('>conitg1\n%s') % (aswat))
  
  elif Nheterozygots > 0 
 
print("Assembly finished with "+str(assembled_N_samam)+" successfully assembled samam from "+str(total_N_samam))

   
  

   
	
	
