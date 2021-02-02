#!/usr/bin/env python3
import sys, os, re, argparse
from agouti import agouti_alleles
from SINE import agouti_SINE
from extension import extension_alleles
from dominantBlack import dominantBlack_alleles
from dilute import dilute_alleles



###--------- USER INPUT ---------###
#Initiate argument parser
parser = argparse.ArgumentParser() 		

#Add arguments
parser.add_argument("--bam", "-i", help="Path to bamfile") 			
parser.add_argument("--out", "-o", help="Name output file")		

#Read arguments from command line
args = parser.parse_args() 	
bam = args.bam
output = args.out


###--------- INITIALIZE --------###

#Remove previous outfile and open new outfile
if os.path.exists(output+'.txt'):
	os.remove(output+'.txt')
outfile = open(output+'.txt', 'a')


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


###--------- CALL FUNCTIONS ---------###
agouti_alleles(bam, outfile, translationDict)
agouti_SINE(bam, outfile)
extension_alleles(bam, outfile, translationDict)
dominantBlack_alleles(bam, outfile, translationDict)
dilute_alleles(bam, outfile, translationDict)

outfile.close()
