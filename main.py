#!/usr/bin/env python3
import sys, os, re
from agouti import agouti_alleles, agouti_SINE
from extension import extension
from dilution import dilution
from dominantBlack import dominantBlack



#Input sample name
if len(sys.argv) == 1:		# no commandline arguments
    sample = input("Please enter a filename: ")
elif len(sys.argv) == 2:	   # something is there
    sample = sys.argv[1]
else:
    sys.stderr.write("Usage: #filename.py <filename>\n")
    sys.exit(1)




ID = sample.split('.')[0]

bamfile = open('/groups/hologenomics/shyam/data/forJulie/bams/'+sample+'.realigned.bam','r')

os.system('rm colourSummary.'+ID+'.txt') #Delete previous outfile	
outfile = open('colourSummary.'+ID+'.txt', 'a')

translationDict = {'TTT':'F','TTC':'F',
					'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
					'ATT':'I','ATC':'I','ATA':'I',
					'ATG':'M',
					'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
					'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
					'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
					'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
					'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
					'TAT':'Y','TAC':'Y',
					'CAT':'H','CAC':'H',
					'CAA':'Q','CAG':'Q',
					'AAT':'N','AAC':'N',
					'AAA':'K','AAG':'K',
					'GAT':'D','GAC':'D',
					'GAA':'E','GAG':'E',
					'TGT':'C','TGC':'C',
					'TGG':'W',
					'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
					'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
					'TAA':'*','TAG':'*','TGA':'*'} #Stop codon


#Call functions
agouti_alleles(sample, ID, bamfile, outfile, translationDict)
agouti_SINE(sample, ID, bamfile, outfile)
extension(sample, ID, bamfile, outfile, translationDict)
dominantBlack(sample, ID, bamfile, outfile, translationDict)
dilution(sample, ID, bamfile, outfile, translationDict)

outfile.close()