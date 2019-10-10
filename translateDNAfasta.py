#! /usr/bin/env python
# Protein translator
import re
import os
 
def BaseChecker(Line):
    ValidBases = set(['A', 'C', 'G', 'T'])
    BaseList = set(Line)
    if BaseList <= ValidBases: return True else: return False def Translator(Sequence): print "For sequence #%s: " % LineNumber if len(Line) % 3 != 0: print "Warning: your sequence is not divisible by 3! Some nucleotides will remain untranslated." for frame in range(0,3): # Translate sequence into a list of codons CodonList = [ ] for x in range(frame, len(Sequence), 3): CodonList.append(Sequence[x:x+3]) # For each codon, translate it into amino acid, and create a list ProteinSeq = [ ] for codon in CodonList: if codon in CodonDict: ProteinSeq.append(CodonDict[codon]) else: break # For each amino acid, calculate its weight ProteinWeight = [ ] for codon in ProteinSeq: ProteinWeight.append(AminoDict[codon]) # Print result, for a given sequence, for a given frame print "Translated in frame %d: %s (%.1f Da)" % ((frame+1), ''.join(ProteinSeq), sum(ProteinWeight)) # Check position of stop codon, making sure it's at the end, and the only one XCount = 0 XEnd = re.search(r"X$", ''.join(ProteinSeq)) for acid in ProteinSeq: if acid == "X": XCount += 1 if XCount == 1 and XEnd is not None: print " Good job: only one stop codon, and at the end!" elif XCount == 1 and XEnd is None: print " WARNING: Stop codon not at end!" elif XCount > 1:
            print " WARNING: multiple stop codons found!"
        elif XCount == 0:
            print " WARNING: No stop codon found!"
    print ""
 
AminoDict={'A':89.09,   'R':174.20, 'N':132.12, 'D':133.10, 'C':121.15, 
'Q':146.15, 'E':147.13, 'G':75.07,  'H':155.16, 'I':131.17, 'L':131.17, 
'K':146.19, 'M':149.21, 'F':165.19, 'P':115.13, 'S':105.09, 'T':119.12, 
'W':204.23, 'Y':181.19, 'V':117.15, 'X':0.0,    '-':0.0,    '*':0.0}
 
CodonDict={'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',  
'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',  
'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',  
'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',  
'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',  
'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',  
'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',  
'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',  
'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',  
'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',  
'AGA':'R',  'AGG':'R',  'TAA':'X',  'TAG':'X',  'TGA':'X'}

fastalist=[]
for f in os.listdir("."):
    if f.endswith(".fasta"):
        fastalist.append(f)
for p in fastalist:
    DNASeq = open(p, 'r')
    LineNumber = 1
    for Line in DNASeq:
        Line = Line.strip('\n').upper().replace(" ","")
        if BaseChecker(Line) == True:
            Translator(Line)    # This is the command that translates sequence by sequence
            LineNumber += 1
        else:
            print "For sequence #%s:\nInvalid bases detected!\n" % LineNumber
            LineNumber += 1
    DNASeq.close()
