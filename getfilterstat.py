#!/usr/bin/python
#in the folder with stderr from filtering job P2 run 
#grep Input *e135* > filterstat.txt (could implement it also in this script)
#this script will transform the txt file into tsv
import itertools
with open("filterstat.tsv","w") as f:
	f.write('Name\tInput_pairs\tboth_survived\tbs_perc\tF_survived\tFs_perc\tR_survived\tRs_perc\tdropped\tdrop_perc\n')
	file=open("trimlog.txt","r")
	lines = file.readlines()
	stalines = lines[10::13]
	names = lines[12::13]
	#print (stalines, names)
	for name,l in zip(names,stalines):
		s = l.split(" ")
		name = name.strip('\n')
		Input_pairs, both_surv, bs_perc, F_surv, Fs_perc, R_surv, Rs_perc, drop, drop_perc = s[3],s[6],s[7],s[11],s[12],s[16],s[17],s[19],s[20]
		#print (name)
		#print (l)
		#print(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n') % (name, Input_pairs, both_surv, bs_perc, F_surv, Fs_perc, R_surv, Rs_perc, drop, drop_perc))
		f.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s') % (name, Input_pairs, both_surv, bs_perc, F_surv, Fs_perc, R_surv, Rs_perc, drop, drop_perc))
	file.close()
