#!/usr/bin/env python3
import sys, os, re


def agouti_SINE(bam, outfile):
	outfile.write('#------------------------SINE insertion---------------------------#\n')
	
	#Initialize files
	aFragmentFile = (os.path.join(os.path.dirname(__file__), 'aFragmentReference.fa'))
	bFragmentFile = (os.path.join(os.path.dirname(__file__), 'bFragmentReference.fa'))	
	intron1A_bed = (os.path.join(os.path.dirname(__file__), 'intron1A.bed'))

	
	#Make intron1A consensus sequence
	os.system('angsd -i '+bam+' -doCounts 1 -doFasta 2 -r chr24: -out chr24.consensus')
	os.system('gunzip chr24.consensus.fa.gz')
	os.system('bedtools getfasta -fi chr24.consensus.fa -fo intron1A.consensus.fa -bed '+intron1A_bed) 
	os.system('rm chr24.consensus.fa chr24.consensus.arg chr24.consensus.fa.fai')
	

	###------------ fragment a ------------###
	#Local alignment
	os.system('water intron1A.consensus.fa '+aFragmentFile+' -gapopen 10.0 -gapextend 0.5 -outfile aFragment.water')
	aFragmentFile = open('aFragment.water', 'r')
	
	#Find alignment positions and stats
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
	
	
	outfile.write('\n')
	

	###------------ fragment b ------------###
	#Local alignment
	os.system('water intron1A.consensus.fa '+bFragmentFile+' -gapopen 10.0 -gapextend 0.5 -outfile bFragment.water')
	bFragmentFile = open('bFragment.water', 'r')
	
	#Find alignment positions and stats
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
	
	
	#Clean-up
	os.system('rm aFragment.water')
	os.system('rm bFragment.water')
	os.system('rm intron1A.consensus.fa')
	
	
	outfile.write('\n')
	
	
	#Distance between a and b fragments
	distance = int(bFragmentStart) - int(aFragmentEnd)
	
	#Print output
	outfile.write('a_start\tb_end\tdistance\n')
	outfile.write(aFragmentEnd+'\t'+bFragmentStart+'\t'+str(distance)+'\n')
	outfile.write('\n')
	
	print('> SINE done')
