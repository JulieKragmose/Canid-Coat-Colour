#!/usr/bin/env python3
import sys, os

def intensity(sample, ID, bamfile, outfile, translationDict):
	
	#Print alignment for all samples
	"""
	samplefile = open('/groups/hologenomics/juliej/data/marta/MHC/references/sampleNames.txt','r')
	samples = []
	for line in samplefile:
		sample = line.rstrip()
		samples.append(sample)
	samplefile.close()
	
	for sample in samples:
		print(sample)
		os.system('samtools tview -d T -p chr20:55855085 -s STR /groups/hologenomics/shyam/data/forJulie/bams/'+sample+'.realigned.bam')
		print('\n')
	"""
	
	outfile.write('###---------------------Intensity (MFSD12)----------------------###\n')
	
	#Initialize
	allelePositions = {'51': '55855085'}
					
	for AApos, BPpos in allelePositions.items():
		outfile.write('-------\n')
		outfile.write('AA '+str(AApos)+'\n')
		outfile.write('-------\n')
		
		
		#Make file with aligment at allele position
		try:
			tview = 'samtools tview -d T -p chr20:'+BPpos+' -s STR /groups/hologenomics/shyam/data/forJulie/bams/'+sample+'.realigned.bam > MFSD12.'+ID+'.pos'+AApos+'.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
		
		
		#Go through codon for each read
		viewfile = open('MFSD12.'+ID+'.pos'+AApos+'.txt', 'r')
		codons = []
		lineNo = 0
		for line in viewfile:
			lineNo += 1
			if lineNo > 3:
				codon = line[:3].strip()
				codons.append(codon)
				if codon != '':
					print(line[:3])
		viewfile.close()
		os.system('rm MFSD12.*') #Clean-up
		
		
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
	
	
	outfile.write('\n')
	print('> Intensity done')