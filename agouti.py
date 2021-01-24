#!/usr/bin/env python3
import sys, os, re


def agouti_alleles(sample, ID, bamfile, outfile, translationDict):
	outfile.write('#-------------------------Agouti (ASIP)---------------------------#\n')
	
	allelePositions = {'82': '23393510', '83': '23393513', '96': '23393552'}
	

	#check codon for each allele position
	for AApos, BPpos in allelePositions.items():
		
		outfile.write('-------\n')
		outfile.write('AA '+str(AApos)+'\n')
		outfile.write('-------\n')
		
		
		#Make file showing aligment at allele position
		try:
			tview = 'samtools tview -d T -p chr24:'+BPpos+' -s STR /groups/hologenomics/shyam/data/forJulie/bams/'+sample+'.realigned.bam > ASIP.'+ID+'.pos'+AApos+'.txt'
			if os.system(tview) != 0:
				raise Exception()
		except:
			sys.exit(1)
		
		
		#Go through codon for each read
		viewfile = open('ASIP.'+ID+'.pos'+AApos+'.txt', 'r')
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
		os.system('rm ASIP.*') #Clean-up
		
		#Count unique codons at allele position
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
		
		if readCount == 0:						#No coverage
			outfile.write('No coverage!\n')
			
		else:									#Coverage
			outfile.write('\n')
			for codon, count in codonDict.items():
				if codon in translationDict:
					AA = translationDict[codon]	#Translate into amino acid
					outfile.write(AA+'\t'+codon+'\t'+str(count)+'/'+str(readCount)+'\n')
			outfile.write('\n')

	print('> Agouti done')



#-------------------------------SINE insertion---------------------------------#
def agouti_SINE(sample, ID, bamfile, outfile):
	outfile.write('#------------------------SINE insertion---------------------------#\n')
	
	
	#Make intron1A consensus sequence
	os.system('angsd -i /groups/hologenomics/shyam/data/forJulie/bams/'+sample+'.realigned.bam -doCounts 1 -doFasta 2 -r chr24: -out chr24.'+ID+'.consensus')
	os.system('gunzip chr24.'+ID+'.consensus.fa.gz')
	bedfile = '/groups/hologenomics/juliej/data/marta/scripts/colour_program/intron1A.bed'
	os.system('bedtools getfasta -fi chr24.'+ID+'.consensus.fa -fo intron1A.'+ID+'.consensus.fa -bed '+bedfile) 
	os.system('rm chr24.'+ID+'.*')
	aFragmentFile = '/groups/hologenomics/juliej/data/marta/scripts/colour_program/aFragmentReference.fa'
	bFragmentFile = '/groups/hologenomics/juliej/data/marta/scripts/colour_program/bFragmentReference.fa'
	
	
	#Local alignment of fragment a
	os.system('water intron1A.'+ID+'.consensus.fa '+aFragmentFile+' -gapopen 10.0 -gapextend 0.5 -outfile aFragment.'+ID+'.water')
	aFragmentFile = open('aFragment.'+ID+'.water', 'r')
	
	#Inititalize
	aPositions = []
	flag = False
	endPositionRegex = re.compile(r'[TAGCtagcN-]+\s+(\d+)\n')
	
	outfile.write('a fragment local alignment:\n')
	for line in aFragmentFile:
		
		#Print stats
		if line.startswith('# Length:'):
			flag = True
		if flag is True:
			outfile.write(line[:-1]+'\n')
		if line.startswith('# Score:'):
			flag = False
			
		#Find end position of a fragment
		if not line.startswith('#') and not line.startswith(' '):
			searchResult = endPositionRegex.search(line)
			if searchResult is not None:
				aPositions.append(searchResult.group(1))
	aFragmentFile.close()
	aFragmentEnd = aPositions[-2]
	os.system('rm aFragment.'+ID+'.water')
	
	
	outfile.write('\n')
	
	
	#Local alignment of fragment b
	os.system('water intron1A.'+ID+'.consensus.fa '+bFragmentFile+' -gapopen 10.0 -gapextend 0.5 -outfile bFragment.'+ID+'.water')
	bFragmentFile = open('bFragment.'+ID+'.water', 'r')
	
	#Initialize
	bPositions = []
	flag = False
	startPositionRegex = re.compile(r'(\d+)\s{1}[TtGgCcAaN-]+')
	
	outfile.write('b fragment local alignment:\n')
	for line in bFragmentFile:
		
		#Print stats
		if line.startswith('# Length:'):
			flag = True
		if flag is True:
			outfile.write(line[:-1]+'\n')
		if line.startswith('# Score:'):
			flag = False
		
		#Find start position of b fragment
		if not line.startswith('#') and not line.startswith(' '):
			searchResult = startPositionRegex.search(line)
			if searchResult is not None:
				bPositions.append(searchResult.group(1))
	bFragmentFile.close()
	bFragmentStart = bPositions[0]
	os.system('rm bFragment.'+ID+'.water')
		
	os.system('rm intron1A.'+ID+'.consensus.fa') #Clean-up
	
	
	outfile.write('\n')
	
	
	#Distance between a and b fragments
	distance = int(bFragmentStart) - int(aFragmentEnd)
	
	#Print output
	outfile.write('a_start\tb_end\tdistance\n')
	outfile.write(aFragmentEnd+'\t'+bFragmentStart+'\t'+str(distance)+'\n')
	outfile.write('\n')
	
	print('> SINE done')
