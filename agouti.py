#!/usr/bin/env python3
import sys, os


def agouti_alleles(bam, outfile, translationDict):
	outfile.write('#-------------------------Agouti (ASIP)---------------------------#\n')


	allelePositions = {'82': '23393510', 
					   '83': '23393513', 
					   '96': '23393552'}


	#Go through each allele
	for AApos, BPpos in allelePositions.items():
		
		
		outfile.write('-------\n')
		outfile.write('AA '+str(AApos)+'\n')
		outfile.write('-------\n')
		
		
		#Make file with aligment at allele position
		try:
			tview = 'samtools tview -d T -p chr24:'+BPpos+' '+bam+' > ASIP.pos'+AApos+'.tview.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
		
		
		#Go through codon for each read
		viewfile = open('ASIP.pos'+AApos+'.tview.txt', 'r')
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
		os.system('rm ASIP.pos'+AApos+'.tview.txt') #Clean-up
		
		#Count unique codons at allele position
		codonDict = dict()
		for codon in codons:
			codon = (codon.strip()).upper()
			if codon != '':
				if codon in codonDict:
					codonDict[codon] += 1
				else:
					codonDict[codon] = 1
					
		#Check coverage
		readCount = 0
		for count in codonDict.values():
			readCount += count
		
		#No coverage
		if readCount == 0:						
			outfile.write('No coverage!\n')
		#There is coverage
		else:					
			outfile.write('\n')
			for codon, count in codonDict.items():
				if codon in translationDict:
					AA = translationDict[codon]	#Translate into amino acid
					outfile.write(AA+'\t'+codon+'\t'+str(count)+'/'+str(readCount)+'\n')
			outfile.write('\n')

	print('> Agouti done')