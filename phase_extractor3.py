#! /usr/bin/env python
# this script parses output file of PHASE
# it creates table with summary statistics per populations, fasta file per population, nexus file with all sequences and haplotypes only (can be loaded to popart) and trait file for popart
# supply .out file from PHASE as the first argument of this script
# written by Jakub Vlcek may 2017


import sys # for using arguments
import re
from collections import Counter
import csv
import dendropy # module for haplotype stats
input_file=sys.argv[1]
#print(input_file)

inp=open(input_file, "rU") # load input file
lines=inp.readlines() # make list of lines
#print(lines)
#print(inp)
def getline(linestr):
	e=list(filter(lambda x:linestr in x,enumerate(lines)))[0][0]  # this is better way how to get line number where specific string occurs and make func out of it
	return e
lineiter=0
lhb=getline("BEGIN BESTPAIRS_SUMMARY\n")
lhe=getline("END BESTPAIRS_SUMMARY\n")

#print(lines[47])
locusname = input_file.split("_")[0] # locus name from the input file, first part separated by "_"
pops=[] #empty lists
data=[]
allist=[]
for lh in lines[lhb+1:lhe]: # this loop generates a set of populations and all data on genotype (only numbers not sequences)
	#print(lh)
	idname=lh.split("(")[0] # split line with samplename and population name
	haps=lh.split("(")[1].strip("\n") # get the haplotypes
	pop=idname.split("_")[2].strip(" ") # get population only
	idn=idname.split("_")[1] # get samplename only
	pop=pop.replace(":","")
	haps=haps.replace(")","")
	pops.append(pop)
	line=[idn,pop,haps] # just make list from the line of genotype data so I do not need to repeat it again	
	allist.append(line) # append every line to allist
	
# haplotype phasing quality (check how many possitions has been resolved poorly and print it to stdout)

phaseprobb = getline("BEGIN PHASEPROBS\n") 
phaseprobe = getline("END PHASEPROBS\n")
lowprob = {}
for ln, lh in enumerate (lines[phaseprobb+1:phaseprobe]):
	probs = re.findall("\d.\d\d",lh) # get the probabilities from the line
	prob = [float(x) for x in probs if float(x) < 0.8] # select the probs that are lower than 0.8	
	if prob != []:	
		lowprob[allist[ln][0]]=prob # append the low probs to the list, maybe I could exclude those samples from subsequent analysis?
print ("positions with phase probability lower than 0.8:")
print (lowprob)

# sequence summary generator
seqlineb=getline("BEGIN LIST_SUMMARY\n")
seqlinee=getline("END LIST_SUMMARY\n")
seq_dict={}

for lh in lines[seqlineb+1:seqlinee]:
	#print (lh)
	hap_id=re.split("\s+",lh)[1]
	hap_seq=re.split("\s+",lh)[2]
	seq_dict[hap_id]=hap_seq # dictionary haplotype names as keys and sequences as values.
	
setpop=list(set(pops)) # list unique populations only (set does the trick)
#print((setpop))
popdict={}
sethap=set()
pop_summary={}# dictionary containing population summary
pop_seq_dict={}
taxlabelsa=[] # lists for nexus file writing
taxlabelsb=[]
nexseqa=[]
nexseqb=[]

for pn, p in enumerate(setpop): # for each individual population
	pophaps=[]
	infolist=[] # list of information for each population (this will form a row in pop_summary)
	indnum=0
	hh=0
	pop_seq=[]
	fasta_name=p+".fasta"
	fasta_pop=open(fasta_name,"w") # create a file with unique haplotypes per population
	for lnn, lh in enumerate(allist): # go through each line (each individual)
		if lh[1] == p:
			#print(lh)
			indnum+=1
			hap1=lh[2].split(",")[0]
			hap2=lh[2].split(",")[1]
			sethap.add(hap1) # generate haplotype set (unique haplotypes per population)
			sethap.add(hap2) # it take both haplotypes
			
			taxlabelsa.append("\t%s_A%s\n" % (lh[0],hap1)) # for nexus taxa names
			taxlabelsb.append("\t%s_B%s\n" % (lh[0],hap2))
			nexseqa.append("\t%s_A%s %s\n" % (lh[0],hap1,seq_dict[hap1])) # for nexus actual sequences
			nexseqb.append("\t%s_B%s %s\n" % (lh[0],hap2,seq_dict[hap2]))

			fasta_pop.write(">%s_A%s\n%s\n" % (lh[0],hap1,seq_dict[hap1]))#create a fasta file per population that contains all copies of each allele (2 copies per individual) this file will be later analysed by dendropy
			fasta_pop.write(">%s_B%s\n%s\n" % (lh[0],hap2,seq_dict[hap2]))
			pophaps.append(hap1)
			pophaps.append(hap2)
			if hap1 != hap2:
				hh+=1 #number of heterozygots			
				#print(haps)
			else:
				hh+=0
	Hobsperpop=float(hh)/indnum # obs heterozygosity 1.91 (frequency of heterozygots in the population)
	Hobsperpop=str("%.2f"% Hobsperpop)	
	infolist.append(p)	
	infolist.append(Hobsperpop)
	infolist.append(indnum)
	pop_summary[p]=infolist # put into the summary dict
	# number of individuals in the summary dict			
	#print(pophaps)
	popdict[p]=pophaps # dictionary of populations as a key and list of all haps as a corresponding value
	hp=list(set(pophaps))# unique haplotypes per given population
	#print(seq_dict[hap1])
	
	fasta_pop.close()
	


        # analyse sequencial diversity in each population, changed
		
	fasta_name=p+".fasta"
	seqs=dendropy.DnaCharacterMatrix.get(
		path=fasta_name,
		schema="fasta")
	if (len(hp)) > 1:
		nucdiv=dendropy.calculate.popgenstat.nucleotide_diversity(seqs)
		K=dendropy.calculate.popgenstat.num_segregating_sites(seqs)
		TajD=dendropy.calculate.popgenstat.tajimas_d(seqs) #new
		Theta=dendropy.calculate.popgenstat.wattersons_theta(seqs)#new		
		#print(nucdiv,K)
		nucdiv=str("%.4f"% nucdiv)
		TajD=str("%.4f"% TajD)
		Theta=str("%.4f"% Theta)
		pop_summary[p].append(nucdiv)
		pop_summary[p].append(K)
		pop_summary[p].append(TajD)
		pop_summary[p].append(Theta)
	else:
		pop_summary[p].append("NA")
		pop_summary[p].append("NA")
		pop_summary[p].append("NA")
		pop_summary[p].append("NA")
#print(pop_seq_dict['Pinta'])
#print(taxlabelsa)
seqlen = (len(seq_dict['1']))
# nexus file for all copies of alleles
with open("allpop"+locusname+".nexus", "w") as nex:
	nex.write("#NEXUS\n")
	nex.write("BEGIN taxa;\n")
	nex.write("\tdimensions ntax=%d;\n" % (2*len(allist)))
	nex.write("\ttaxlabels\n")
	for n, t in enumerate(taxlabelsa):
		nex.write(t)
		nex.write(taxlabelsb[n])	
	nex.write(";\nend;\n\n")
	nex.write("BEGIN characters;\n")
	nex.write("\tdimensions nchar=%d;\n" % seqlen)
	nex.write("\tformat datatype=dna missing=? gap=-;\n")
	nex.write("\tmatrix\n")
	for n, t in enumerate(nexseqa):
		nex.write(t)
		nex.write(nexseqb[n])
	nex.write("\n;\nend;\nBEGIN SETS;\n")
	hapint = 1
	for p in setpop:
		nhap = len(popdict[p])
		hapint = hapint+nhap
		hup = hapint-1	
		#print(p, nhap)
		nex.write("TaxSet %s = %d-%d;\n" % (p,hapint-nhap,hup))
	
	nex.write("end;\nBEGIN CODONS;\n")

	nex.write("\tCODONPOSSET * codonpositions =\n") # change this according reading frame of each locus
	nex.write("\t\t1: 3-"+str(seqlen-2)+"\\3,\n")
	nex.write("\t\t2: 4-"+str(seqlen-1)+"\\3,\n")
	nex.write("\t\t3: 5-"+str(seqlen)+"\\3;\n")
	nex.write("CODESET * UNTITLED = Universal: all;\nend;")
	
	 
NOhap=len(sethap)
NOpop=len(setpop)
aldist=[[0 for j in range(NOpop)] for i in range(NOhap)] # just make an empty array with columns corresponding to number of pops and rows corresponding to number of haps 

# statistics for all the dataset
allseqs=dendropy.DnaCharacterMatrix.get(path="allpop"+locusname+".nexus", schema="nexus")
nucdivall=dendropy.calculate.popgenstat.nucleotide_diversity(allseqs)
Kall=dendropy.calculate.popgenstat.num_segregating_sites(allseqs)
TajDall=dendropy.calculate.popgenstat.tajimas_d(allseqs) #new
Thetaall=dendropy.calculate.popgenstat.wattersons_theta(allseqs)#new		
print("Total: Pi:%0.6f\tK:%d\tTajD:%0.6f\tTheta/sequence:%1.6f" % (nucdivall,Kall,TajDall,Thetaall))

			
for np, p in enumerate(setpop): # very cool looping trick how to get index of item and the item itself np=index of population and p=population name
	#print(p,np)
	hsum=0
	z=(popdict[p]) # list of haplotypes for given population
	num=Counter(z)	# occurence counter, number of each haplotype in list z
	#print(num)
	hp_div_pop=len(num) # number of haplotypes per population
	pop_summary[p].append(hp_div_pop) # add it to the pop_summary_dict number of haplotypes
	pop_summary[p].append('c'+(str(list(set(z)))))# list of haplotypes 
	for h in range(NOhap): # loop for filling the aldist array
		hs=h+1 # because it starts with 0 but haplotype names start from 1		
		hstr=str(hs) 	 	
		s=num[hstr] # the counter format approached the same as dictionary, s=number of given haplotype h in given population p, if the haplotype is not present it generates 0
		
		hnn=float(s)/(2*pop_summary[p][2]) # realtive allele frequency in the population the (2*pop_summary[p][2]) represents number of individuals * 2 = no chromosomal copies, Nei's gene diversity 
		d=hnn**2 #nove
		hsum+=d
		#print(s)
		#print(np,h,hs,s)
		#nh=int(h)
		aldist[h][np]=s # fill the array
		#outstring="s%\t
	hsexp=1-hsum # NEi He_exp
	hsexp=str("%.2f"% hsexp)
	#print(p,hsum,hsexp)
	pop_summary[p].append(hsexp)	
#print(popdict)

	



with open ("popart_traits"+locusname+".txt", "wb") as f: # export the array by cool module csv stright into csv file
	writer=csv.writer(f, delimiter="\t")
	writer.writerow(setpop) # add header row population names, for popart it needs only to add the names of haplotypes eg in excel 
	writer.writerows(aldist)

	
with open("summary_pop"+locusname+".txt", "w") as i:
	i.write("pop,Hs_obs,Noind,Pi,K,TajD,Theta,Ap,haps,Hs_exp\n")	
	for p in setpop:
		line=str(pop_summary[p])
		line=re.sub("^\[","",line)
		line=re.sub("\]$","",line)
		line=line.replace("[","(").replace("]",")")
		i.write(line+"\n")

with open("haplolist"+locusname+".nexus", "w") as hl: # write nexus file of haplotypes for popart
	hl.write("#NEXUS\n")
	hl.write("begin taxa;\n")
	hl.write("\tdimensions ntax=%d;\n" % NOhap)
	hl.write("\ttaxlabels\n")
	for hap in range(NOhap):
		hap+=1
		hl.write("\t%d\n" % hap)
	hl.write(";\nend;\n\n")
	hl.write("begin characters;\n")
	seqlen=len(seq_dict['1'])
	hl.write("\tdimensions nchar=%d;\n" % seqlen)
	hl.write("\tformat datatype=dna missing=? gap=-;\n")
	hl.write("\tmatrix\n")
	for hap in range(NOhap):
		hap+=1
		hap=str(hap)
		seq=seq_dict[hap]
		hl.write("\t%s %s\n" % (hap,seq))
	hl.write("\n;\nend;")
inp.close()

