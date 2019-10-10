#! /usr/bin/env python

import re


seqfile=open('pro_phase_TLR1B.fasta','r')
popfile=open('finalpopTLR1B.tsv','r')
#num_lines = sum(1 for line in infile)
#print num_lines
LN=0
nameline=range (1,638,3) #num_lines - 1 # read it in text file
firstline=range (2,639,3) #num_lines
secondline=range (3,640,3) #num_lines + 1

popind={} # dictionary
for ind in popfile:
	#print ind			
	indname, amplname = ind.split(" ")[0], ind.split(" ")[1]
	pop=ind.split(" ")[2].strip("\n")
	
	pi = indname+'_'+pop
	popind[amplname]=pi
#print(popind)


for l in seqfile:
	l=l.strip('\n')
	LN=LN+1
	if LN in nameline:
		#print LN
		seqnames=l
		p=popind[seqnames]
		#print "sss: %s" % (seqnames)
		print "%s_%s" % (seqnames, p)
						
		
			
		
			

	elif LN in firstline:
		fl=l.replace('Y','C')
		fl=fl.replace('R','G')
		fl=fl.replace('M','A')
		fl=fl.replace('S','C')
		fl=fl.replace('W','A')
		fl=fl.replace('K','G')
		print fl
	elif LN in secondline:
		sl=l.replace('Y','T')
		sl=sl.replace('R','A')
		sl=sl.replace('M','C')
		sl=sl.replace('S','G')
		sl=sl.replace('W','T')
		sl=sl.replace('K','T')
		print sl
	
nloc=len(fl)
nind=LN/3
locseq=range (1,nloc+1) # it need to replace "[]," by nothing and write P to this line
SNP='S' * nloc
print nind
print nloc
print locseq
print SNP 
seqfile.close()
popfile.close()
