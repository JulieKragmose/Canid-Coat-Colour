#!/usr/bin/env python3
import sys, os

def dominantBlack_alleles(bam, outfile, translationDict):
	outfile.write('###------------------Dominant Black (CBD103)--------------------###\n')
	outfile.write('NB! CBD103 is on the reverse strand!\n')
	
	#Initialize
	allelePositions = {'23': '58965448'}
	
					
	#check codon for each allele position
	for AApos, BPpos in allelePositions.items():
		outfile.write('---------------------\n')
		outfile.write('AA '+str(int(AApos)-2)+' - '+str(int(AApos)+2)+' alignment\n')
		outfile.write('---------------------\n')
		
		
		#Make file with aligment at allele position
		BPpos = int(BPpos) - 6
		try:
			tview = 'samtools tview -d T -p chr16:'+str(BPpos)+' '+bam+' > CBD103.pos'+AApos+'.tview.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
		
		
		#Go through reads for AA 21-25 and save codons
		viewfile = open('CBD103.pos'+AApos+'.tview.txt', 'r')
		indexDict = {'15':[], '12':[], '9':[], '6':[], '3':[],}
		lineNo = 0
		
		for line in viewfile:
			lineNo += 1
			if lineNo > 3:
				if line[:15].strip() != '':
					outfile.write(line[:15]+'\n')
				
				for i in range(15,0,-3):
					reverseCodon = line[i-3:i].strip()
					if reverseCodon != '':
						indexDict[str(i)].append(reverseCodon)
		viewfile.close()
		os.system('rm CBD103.pos'+AApos+'.tview.txt') #Clean-up

		outfile.write('\n')
		
		#Reverse complement the codons
		position = 21
		for reverseCodons in indexDict.values():
			codonDict = dict()
			
			for c in reverseCodons:
				#Reverse
				forwardCodon = c[::-1].upper()
				
				#Complement
				codon = ''
				if forwardCodon in translationDict: #Check if the codon is valid
					for base in forwardCodon:
						if base == 'A':
							codon += 'T'
						elif base == 'T':
							codon += 'A'
						elif base == 'G':
							codon += 'C'
						elif base == 'C':
							codon += 'G'
						else:
							print('ERROR: Invalid nucleotide: '+base)
							sys.exit(1)	
				
					#Count unique codons
					if codon != '':
						if codon in codonDict:
							codonDict[codon] += 1
						else:
							codonDict[codon] = 1
			
			
						
			#Check coverage and print output
			outfile.write('-------\n')
			outfile.write('AA '+str(position)+'\n')
			outfile.write('-------\n')
			
			if len(codonDict) != 0: #Coverage
				for codon, count in codonDict.items():
					if codon in translationDict:
						AA = translationDict[codon]	#Translate into amino acid
						outfile.write(AA+'\t'+codon+'\t'+str(count)+'/'+str(len(reverseCodons))+'\n')
				outfile.write('\n')
				
			else: #No coverage				
				outfile.write('No coverage!\n')
				outfile.write('\n')		
			
			#Next codon
			position += 1
	
	outfile.write('\n')
	
	print('> Dominant black done')