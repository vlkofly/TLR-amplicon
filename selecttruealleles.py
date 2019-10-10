#! /usr/bin/env python
# script that takes collapsed amplicons per individual
# and selects true variants according read frequency
# written by Jakub Vlcek may 2017

import os
import itertools

def getfreq(listofnames): # get the frequency from sequence name
 freq = []
#order = []
 for n in listofnames:
  f = n.split('_')[3].strip('\n') # collect frequency
  #o = n.split('_')[1].strip('\n')[1] # collect order
  freq.append(f)
  #order.append(o)
 return(freq)
 

# collect files that end with fastq in local directory
infasta = []
for f in os.listdir("."):
 if f.endswith('.fastq'):
  infasta.append(f)
infasta.sort()

# loop through all the files and process the frequency
for f in infasta:
 ft = f.replace("collapsed","").replace("assembled","")
 print("processing"+ft)
 with open('../true/true'+ft,'w') as outfile:
  with open(f,'r') as fi:
   lines = fi.readlines()
   seqname = lines[0:8:2]
   seq = lines[1:9:2]
   freq = getfreq(seqname)
   #print(freq)
   #print (float(freq[1])/float(freq[0]))
   if (float(freq[1])/float(freq[0])) > 1/3.0: # if the ratio of second sequence to first is higher than 1/3 we consider the individual heterozygote in this amplicon
    print ("heterozygot")
    outfile.write(lines[0:3])
    outfile.write (('%s\n%s\n%s\n%s') % (seqname[0].strip('\n'), seq[0].strip("\n"), seqname[1].strip('\n'), seq[1].strip("\n")))
    outfile.close 
   else:
    print ("homozygot \n")    
    outfile.write (('%s\n%s') % (seqname[0].strip('\n'), seq[0].strip("\n")))
    outfile.close  
   
	
	
