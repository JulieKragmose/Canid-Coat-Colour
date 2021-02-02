#!/usr/bin/env python3
import sys, os

def extension_alleles(bam, outfile, translationDict):	
	outfile.write('###---------------------Extension (MC1R)------------------------###\n')
	outfile.write('NB! MC1R is on the reverse strand!\n')
	
	#Initialize
	allelePositions = {'264': '63694458', 
					   '78': '63695016', 
					   '306': '63694332', 
					   '301': '63694347'}
					
	#check codon for each allele position
	for AApos, BPpos in allelePositions.items():
		outfile.write('-------\n')
		outfile.write('AA '+str(AApos)+'\n')
		outfile.write('-------\n')
		
		
		#Make file with aligment at allele position
		try:
			tview = 'samtools tview -d T -p chr5:'+BPpos+' '+bam+' > MC1R.pos'+AApos+'.tview.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
		
		
		#Go through each read in alignment
		viewfile = open('MC1R.pos'+AApos+'.tview.txt', 'r')
		codons = []
		lineNo = 0
		for line in viewfile:
			lineNo += 1
			if lineNo > 3:
				reverseCodon = line[:3].strip()
				
				#Print alignment before reverse complementing
				if reverseCodon != '':
					outfile.write(line[:3]+'\n')
					
				#Reverse codon
				forwardCodon = reverseCodon[::-1].upper()
				
				#Complement codon
				if forwardCodon in translationDict: #Check if the codon is valid
					codon = ''
					for base in forwardCodon:
						if base == 'A':
							codon += 'T'
						elif base == 'T':
							codon += 'A'
						elif base == 'G':
							codon += 'C'
						elif base == 'C':
							codon += 'G'
						elif base == ' ':
							codon += ' '
						else:
							print('ERROR: Invalid nucleotide: '+base)
							sys.exit(1)				
					codons.append(codon) 		
		viewfile.close()
		os.system('rm MC1R.pos'+AApos+'.tview.txt') #Clean-up
						
		
		#Count unique codons
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
			
	print('> Extension done')
