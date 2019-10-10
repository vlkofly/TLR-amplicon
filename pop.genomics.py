#! /usr/bin/env python
# 

import sys # for using arguments
import re
import csv
import dendropy # module for haplotype stats
input_file=sys.argv[1] # input file nexus 
#print(input_file)

# statistics for all the dataset
allseqs=dendropy.DnaCharacterMatrix.get(path=input_file, schema="nexus")
nucdivall=dendropy.calculate.popgenstat.nucleotide_diversity(allseqs)
Kall=dendropy.calculate.popgenstat.num_segregating_sites(allseqs)
TajDall=dendropy.calculate.popgenstat.tajimas_d(allseqs) #new
Thetaall=dendropy.calculate.popgenstat.wattersons_theta(allseqs)#new		
print("Total: Pi:%0.6f\tK:%d\tTajD:%0.6f\tTheta/sequence:%1.6f" % (nucdivall,Kall,TajDall,Thetaall))


