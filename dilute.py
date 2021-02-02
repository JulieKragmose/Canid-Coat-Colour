#!/usr/bin/env python3
import sys, os

def dilute_alleles(bam, outfile, translationDict):
	outfile.write('###----------------------Dilution (MLPH)------------------------###\n')
	
	
	#Initialize
	allelePositions = {'199': '48149423'}
		
	
	for AApos, BPpos in allelePositions.items():
		outfile.write('-------\n')
		outfile.write('AA '+AApos+'\n')
		outfile.write('-------\n')
		
		#Make file with aligment at allele position
		try:
			tview = 'samtools tview -d T -p chr25:'+BPpos+' '+bam+' > MLPH.pos'+AApos+'.tview.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
			
		#Go through codon for each read
		viewfile = open('MLPH.pos'+AApos+'.tview.txt', 'r')
		codons = []
		lineNo = 0
		for line in viewfile:
			lineNo += 1
			if lineNo > 3:
				codon = line[:3].strip()
				codons.append(codon)
				if codon != '':
					outfile.write(line[:3]+'\n')
		viewfile.close()
		os.system('rm MLPH.pos'+AApos+'.tview.txt')
		
		#Count codons
		codonDict = dict()
		for codon in codons:
			codon = (codon.strip()).upper()
			if codon != '':
				if codon in codonDict:
					codonDict[codon] += 1
				else:
					codonDict[codon] = 1
					
		#Check coverage and print output
		readCount = 0
		for count in codonDict.values():
			readCount += count
		
		if readCount == 0:				#No coverage
			outfile.write('No coverage!\n')
			
		else:							#Coverage
			outfile.write('\n')
			for codon, count in codonDict.items():
				if codon in translationDict:
					AA = translationDict[codon]	#Translate into amino acid
					outfile.write(AA+'\t'+codon+'\t'+str(count)+'/'+str(readCount)+'\n')
			outfile.write('\n')
			
	print('> Dilution done')