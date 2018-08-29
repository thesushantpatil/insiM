#!/usr/bin/env python

import argparse
import os
import re
import sys
import random
import itertools
from collections import defaultdict
import operator
from operator import itemgetter
import pysam
import datetime
import time
from time import gmtime, strftime
from array import *

print "############################"
print "insiM initiated.\n"

time_start = datetime.datetime.now().replace(microsecond=0)
print "############################"

def main(assay, bamfile, targetbed, mutfraction, outputfastqname, indellength, insertsequence, mutationtype, genomefasta, amplicon, readlength):
	
	# Checking if parameters and files exist
	if (not assay in list(caseIgnore('amplicon'))) and (not assay in list(caseIgnore('capture'))):
		print " -assay (ASSAY TYPE) should be either 'amplicon' or 'capture' (for Hybrid Cature)"
		sys.exit(1)

	if (not os.path.exists(bamfile)):
		print bamfile + " does not exist!"
		sys.exit(1)
	
	if (not os.path.exists(targetbed)):
		print targetbed + " does not exist!"
		sys.exit(1)
		
	if (not os.path.exists(genomefasta)):
		print genomefasta + " does not exist!"
		sys.exit(1)
	
	if not mutationtype in list(caseIgnore("mix")):
		try:
			mutfract = float(mutfraction)
		except:
			print "Mutation fraction must be a number between 0-1"
			sys.exit(1)
		
		if not 0 <= mutfract <= 1:
			print "Mutation fraction must be a number between 0-1"
			sys.exit(1)
	
	mutType=[]
	for item in ['snv', 'ins', 'del', 'dup', 'mix']:
		mutType.extend(list(caseIgnore(item)))

	if not mutationtype in mutType:
		print "Mutation type must be either snv (Single Nucleotide Variant), ins (Insertion), del (Deletion), dup (Duplication) or mix (Mixed)."
		sys.exit(1)
	
	if not ((mutationtype in list(caseIgnore('snv'))) or (mutationtype in list(caseIgnore('mix')))):
		print "mutationtype is %s" %mutationtype
		try:
			indellen = int(indellength)
		except:
			print "Indel length specified is not an integer > 0. Default length of 10 will be used"
			indellen = 10

		if not indellen > 0:
			print "Indel length specified is not an integer > 0. Default length of 10 will be used"
			indellen = 10

	if insertsequence != None:
		if not validSequence(insertsequence):
			print "Sequence to be inserted contains characters other than A/C/G/T"
			sys.exit(1)
	
	amplicon_flag=False # Is this an amplicon mutation?
	
	if assay in list(caseIgnore("amplicon")) and amplicon != None:
		if (not os.path.exists(amplicon)):
			print amplicon + " does not exist!"
			sys.exit(1)

		else:
			amplicon_flag=True
			amplicon=readAmpBed(amplicon)

	if assay in list(caseIgnore("amplicon")):
		try:
			readlength = int(readlength)
		except:
			print "Read length must be an integer > 0"
			sys.exit(1)

	if mutationtype in list(caseIgnore('mix')):
		fileCheck2(targetbed)
		mixed_mutation_flag=True
	else:
		fileCheck1(targetbed)
		mixed_mutation_flag=False


	#Reading in genome FASTA file...
	print "Reading in genome FASTA file ..."
	genome_input = open(genomefasta, 'r')

	genome_lists = {}
	current_chrom = "null"
	for line in genome_input:
		if line.startswith(">"):
			current_chrom = line.rstrip().split(">")[1]
			genome_lists[current_chrom] = []
		else:
			genome_lists[current_chrom].append(line.rstrip().upper())
	
	genome_input.close()

	print "Merging chromosomes into strings for final genome dict ..."
	genome = {}
	for chrom in genome_lists:
		genome[chrom]= "N" + "".join(genome_lists[chrom])  #This is done to make positions "1" based, so 1 is the first position

	print "Opening input and output files..."
	infile=pysam.AlignmentFile(bamfile, "rb")
	if not outputfastqname:
		outputfilebase = bamfile.rsplit('/', 1)[1].split('.bam')[0] + '.' + mutationtype
	else:
		outputfilebase = outputfastqname

	R1fastq = open(outputfilebase + "_R1_001.fastq", 'w')
	R2fastq = open(outputfilebase + "_R2_001.fastq", 'w')
	
	vcf	= open(outputfilebase + ".insiM.vcf", 'w')
	
	# Writing the header of VCF
	vcf.write("##fileformat=VCFv4.1\n")
	vcf.write("##fileDate="+strftime("%Y-%m-%d %H:%M:%S", gmtime())+"\n")
	vcf.write("##source="+bamfile+"\n")
	vcf.write("##reference=file:///"+genomefasta+"\n") 
	vcf.write("##contig=\n") 
	vcf.write("##phasing=partial\n")
	vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n") 
	vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
	vcf.write("##FILTER=<ID=q10,Description=\"Quality below 10\">\n")
	vcf.write("##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n")
	vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
	vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
	vcf.write("##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n")
	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmy_sample\n") 


	# default <read>:<is filtered>:<control number>:<sample number>
	name_appendage_R1 = "1:N:0:0"
	name_appendage_R2 = "2:N:0:0"

	# Putting all bed entries into a list
	print "Reading in bed entries..."
	bed_entries = []
	bed_in = open(targetbed, 'r')
	for line in bed_in:
		bed_entries.append(line.rstrip())
	print "\t" + str(len(bed_entries)) + " bed entries recorded."

	# 'mutation_spots' will be filled with positions that are the midpoints of each BED entry.  This is where the mutations will be created.
	print "Determining mutation locations based on bed file entry mid-points..."
	
	mutation_spots = []
	
	for entry in bed_entries:
		result = get_mut_location(entry)
		if mutationtype in list(caseIgnore('mix')):
			mutation_spots.append(result + '\t' + entry.split('\t')[3])
		else: 
			mutation_spots.append(result)

	chrList = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" ]
	
	for i in range(0, 24):
    	
    	# The plan here is to copy every read in the bam file into either R1 or R2 "segments".  AlignedSegments is pysam-ese for an aligned read.
		# These are split into R1 and R2 is so it's easier to output them to the final files as pairs (e.g. into paired fastq files)
		# So all reads are written into these dictionaries and then selectively modified by using the pileup function to find mutation spots/reads.
		# After modification of randomly selected reads, the entire read set can then be sent out to output files.
		print "Working on", chrList[i], "..."

		R1_segments = {}
		R2_segments = {}

		for segment in infile.fetch(chrList[i]):
			name = segment.query_name
			if segment.is_secondary == False:
				if segment.is_read1 == True:
					R1_segments[name] = pysam.AlignedSegment()
					R1_segments[name] = segment
				else:
					R2_segments[name] = pysam.AlignedSegment()
					R2_segments[name] = segment
			
		# Iterating through 'mutation_spots' to mutate specific reads in the read1 and read2 dictionaries..."	

		for mutation_loc in mutation_spots:
			mutation_chrom = mutation_loc.split('\t')[0]
			if mutation_chrom == chrList[i]:
				mutation_start = int(mutation_loc.split('\t')[1])
				mutation_end = mutation_start + 1

				if mutationtype in list(caseIgnore("del")) and amplicon_flag==True and ((mutation_start + indellen) >= checkPosBed(mutation_chrom,mutation_start,amplicon)[2]):
					print "Deletion region for the following locus falls outside the amplicon boundry and thus will not be mutated - %s" %mutation_loc
					continue

				if mixed_mutation_flag == True:
					
					mutationtype = str(mutation_loc.split('\t')[2].split(';')[0])
					
					mutfract = float(mutation_loc.split('\t')[2].split(';')[1])
					
					insertsequence = str(mutation_loc.split('\t')[2].split(';')[2])
					if (not len(insertsequence) > 0) or (not validSequence(insertsequence)):
						insertsequence = None

					try:
						indellen = int(mutation_loc.split('\t')[2].split(';')[3])
						if not indellen > 0:
							print "Indel length for mutation locus (%s) is not an integer > 0. Default length of 10 will be used" %mutation_loc
							indellen = 10
					except:
						print "Indel length for mutation locus (%s) is not an integer > 0. Default length of 10 will be used" %mutation_loc
						indellen = 10

				pile = infile.pileup(mutation_chrom, mutation_start, mutation_end, max_depth = 1000000)		# Generating a mutation-locus specific pileup construct of localized reads
				
				for pileupcolumn in pile:
					pileup_reads={}
					if pileupcolumn.pos == mutation_start:
						count = 0 # For counting total reads at a locus
						count2 = 0 # For counting mutated reads at a loucs
				
						for pr in pileupcolumn.pileups:
							count = count + 1
							if pr.alignment.query_name not in pileup_reads:
								pileup_reads[pr.alignment.query_name]=[pr]
							else:
								pileup_reads[pr.alignment.query_name].append(pr)

						if (mutationtype in list(caseIgnore("ins"))) and (insertsequence == None): 
							insertseq = ''.join(random.choice(['A', 'C', 'G', 'T']) for x in range(indellen))
						elif (((mutationtype in list(caseIgnore("snv"))) and (insertsequence == None)) or ((mutationtype in list(caseIgnore("del"))) and (insertsequence == None)) or ((mutationtype in list(caseIgnore("dup"))) and (insertsequence == None))):
							insertseq = None
						elif insertsequence != None:
							insertseq = str(insertsequence) 
				
						# For each pair or a single read (r1/r2)
						for pileupPair in pileup_reads:
							if random.randint(1,100) <= mutfract*100:
								for pileupread in pileup_reads[pileupPair]:
									if not pileupread.is_del and not pileupread.is_refskip:
											count2 = count2 + 1
											if pileupread.alignment.is_read1 == True:
												original_sequence = R1_segments[pileupread.alignment.query_name].query_sequence
												original_qscores = R1_segments[pileupread.alignment.query_name].query_qualities
											else:
												original_sequence = R2_segments[pileupread.alignment.query_name].query_sequence
												original_qscores = R2_segments[pileupread.alignment.query_name].query_qualities

											if mutationtype in list(caseIgnore("snv")):
												mutated_sequence, insertseq = snv(original_sequence, pileupread.query_position, genome, mutation_chrom, mutation_start, insertseq)
											
											elif mutationtype in list(caseIgnore("ins")) and amplicon_flag==False:
												mutated_sequence = insert(original_sequence, pileupread.query_position, insertseq)

											elif mutationtype in list(caseIgnore("ins")) and amplicon_flag==True:
												mutated_sequence = insertAMP(pileupread.alignment, original_sequence, pileupread.query_position, insertseq, genome, mutation_chrom, mutation_start, amplicon, readlength)
												original_qscores = original_qscores[:pileupread.query_position]+array('B',[original_qscores[pileupread.query_position]]*indellen)+original_qscores[pileupread.query_position:]
												if (pileupread.alignment.is_reverse == False):
													original_qscores = original_qscores[:readlength] + array('B',[original_qscores[len(original_qscores)-1]]*(len(mutated_sequence)-len(original_qscores)))
												else:
													original_qscores = array('B',[original_qscores[0]]*(len(mutated_sequence)-len(original_qscores))) + original_qscores[:readlength]
												
											elif mutationtype in list(caseIgnore("del")) and amplicon_flag==False:
												try:
													mutated_sequence = deletion(infile, pileupread.alignment, genome, mutation_chrom, mutation_start, indellen)
													original_qscores = original_qscores[0:len(mutated_sequence)] 
												except:
													continue

											elif mutationtype in list(caseIgnore("del")) and amplicon_flag==True:
												mutated_sequence = deletionAMP(pileupread.alignment, original_sequence, pileupread.query_position, indellen, genome, mutation_chrom, mutation_start, amplicon, readlength)
												original_qscores = original_qscores[:pileupread.query_position]+original_qscores[(pileupread.query_position+indellen):] #trim qc score to stay within amplicon
												if (pileupread.alignment.is_reverse == False):
													original_qscores = original_qscores[:readlength] + array('B',[original_qscores[len(original_qscores)-1]]*(len(mutated_sequence)-len(original_qscores)))
												else:
													original_qscores = array('B',[original_qscores[0]]*(len(mutated_sequence)-len(original_qscores))) + original_qscores[:readlength]
											
											elif mutationtype in list(caseIgnore("dup")) and amplicon_flag==False:
												mutated_sequence, insertseq = duplication(original_sequence, pileupread.query_position, indellen, genome, mutation_chrom, mutation_start)
											
											elif mutationtype in list(caseIgnore("dup")) and amplicon_flag==True:
												mutated_sequence, insertseq = duplicationAMP(pileupread.alignment, original_sequence, pileupread.query_position, indellen, genome, mutation_chrom, mutation_start, amplicon, readlength)
												original_qscores = original_qscores[:pileupread.query_position]+array('B',[original_qscores[pileupread.query_position]]*indellen)+original_qscores[pileupread.query_position:]
												if (pileupread.alignment.is_reverse == False):
													original_qscores = original_qscores[:readlength] + array('B',[original_qscores[len(original_qscores)-1]]*(len(mutated_sequence)-len(original_qscores)))
												else:
													original_qscores = array('B',[original_qscores[0]]*(len(mutated_sequence)-len(original_qscores))) + original_qscores[:readlength]
											else:
												print "It should be impossible to get in here..."
												sys.exit(1)

											# Adding mutated data to fastq dict
											if pileupread.alignment.is_read1 == True:
												R1_segments[pileupread.alignment.query_name].query_sequence = mutated_sequence 
												R1_segments[pileupread.alignment.query_name].query_qualities = original_qscores
											else:
												R2_segments[pileupread.alignment.query_name].query_sequence = mutated_sequence
												R2_segments[pileupread.alignment.query_name].query_qualities = original_qscores

						# Creating VCF entry
						if  ((mutationtype in list(caseIgnore('snv'))) or ((mixed_mutation_flag==True) and (mutationtype in list(caseIgnore('snv'))))):			
							REF = genome[mutation_chrom][mutation_start]
							vcf.write(mutation_chrom + '\t' + str(mutation_start) + '\t' + '.' + '\t' + REF + '\t' + insertseq + '\t' + '.' + '\t' + '.' + '\t')
							vcf.write('DP=' + str(count) + ';AF=' + str(round(float(count2)/count,3)) + ';REF=' + REF + ';ALT=' + insertseq + ';' + '\t' + 'GT:DP' + '\t' + ':' + '\n')

						elif ((mutationtype in list(caseIgnore('del'))) or ((mixed_mutation_flag==True) and (mutationtype in list(caseIgnore('del'))))):
							REF = genome[mutation_chrom][(mutation_start-1):(mutation_start+indellen)]
							ALT = genome[mutation_chrom][mutation_start-1]
							vcf.write(mutation_chrom + '\t' + str(mutation_start - 1) + '\t' + '.' + '\t' + REF + '\t' + ALT + '\t' + '.' + '\t' + '.' + '\t')
							vcf.write('DP=' + str(count) + ';AF=' + str(round(float(count2)/count,3)) + ';REF=' + REF + ';ALT=' + ALT + ';' + '\t' + 'GT:DP' + '\t' + ':' + '\n')

						elif ((mutationtype in list(caseIgnore('ins'))) or (mutationtype in list(caseIgnore('dup'))) or ((mixed_mutation_flag==True) and ((mutationtype in list(caseIgnore('ins'))) or (mutationtype in list(caseIgnore('dup')))))):
							REF = genome[mutation_chrom][mutation_start]
							ALT = genome[mutation_chrom][mutation_start]+insertseq
							vcf.write(mutation_chrom + '\t' + str(mutation_start) + '\t' + '.' + '\t' + REF + '\t' + ALT + '\t' + '.' + '\t' + '.' + '\t')
							vcf.write('DP=' + str(count) + ';AF=' + str(round(float(count2)/count,3)) + ';REF=' + REF + ';ALT=' + ALT + ';' + '\t' + 'GT:DP' + '\t' + ':' + '\n')

		# Outputing results to R1 and R2 FASTQ files
		for readname in R1_segments:
			if readname in R2_segments:
				if R1_segments[readname].is_reverse == True:
					R1sequence = revcomp(R1_segments[readname].query_sequence)
					R1qscores = reversed_string(quality_score_encoder(R1_segments[readname].query_qualities))
				else:
					R1sequence = R1_segments[readname].query_sequence
					R1qscores = quality_score_encoder(R1_segments[readname].query_qualities)
				if amplicon_flag == True:
					R1sequence = R1sequence + "N"*(readlength - len(R1sequence))
					R1qscores = R1qscores + R1qscores[len(R1qscores)-1]*(readlength - len(R1qscores))
				R1fastq.write("@" + readname + ' ' + name_appendage_R1 + '\n')
				R1fastq.write(R1sequence + '\n')
				R1fastq.write("+\n")
				R1fastq.write(R1qscores + '\n')
			
				if R2_segments[readname].is_reverse == True:
					R2sequence = revcomp(R2_segments[readname].query_sequence)
					R2qscores = reversed_string(quality_score_encoder(R2_segments[readname].query_qualities))
				else:
					R2sequence = R2_segments[readname].query_sequence
					R2qscores = quality_score_encoder(R2_segments[readname].query_qualities)
				if amplicon_flag == True:
					R2sequence = R2sequence + "N"*(readlength - len(R2sequence))
					R2qscores = R2qscores + R2qscores[len(R2qscores)-1]*(readlength - len(R2qscores))
				R2fastq.write("@" + readname + ' ' + name_appendage_R2 + '\n')
				R2fastq.write(R2sequence + '\n')
				R2fastq.write("+\n")
				R2fastq.write(R2qscores + '\n')

	vcf.close()
	R1fastq.close()
	R2fastq.close()
	infile.close()

	time_end = datetime.datetime.now().replace(microsecond=0)
	print "insiM completed in %s seconds" %(str((time_end-time_start).total_seconds()))


####################################

############ FUNCTIONS #############

####################################

def configParser(file):
	print "Parsing the configuration file ..."

	try:
		# Read parameters from configuration file in to a dictionary 
		params = {}
		with open(file) as f:
			for line in f:
				line = line.strip()
				if not line.startswith("#"):
					key_value = line.split("=")
					if len(key_value) == 2:
						params[key_value[0].strip()] = key_value[1].strip()

		# Convert the read parameters in to desired types
		params["assay"]	= str(params["assay"])
		params["bam"]	= str(params["bam"])
		params["target"]	= str(params["target"])
		params["mutation"]	= str(params["mutation"])
		params["genome"]	= str(params["genome"])
		if params["ampliconbed"]:
			params["ampliconbed"]	= str(params["ampliconbed"])
		if params["read"]:
			params["read"]	= int(params["read"])
		if params["vaf"]:
			params["vaf"]	= float(params["vaf"])
		if params["len"]:
			params["len"]	= int(params["len"])
		if params["seq"]:
			params["seq"]	= str(params["seq"])
		if params["out"]:
			params["out"]	= str(params["out"])

		return params

	except:
		print "ERROR: Configuration file parsing error"


def fileCheck1(file):
	print "Parsing the target bed file for format check ..."

	with open(file) as f:
		for line in f:
			try:
				start = int(line.split('\t')[1])
				end = int(line.split('\t')[2])
			except:
				print "ERROR: Chromosome postions in target bed file are not integers"
				print line
				print int(line.split('\t')[1])
				print int(line.split('\t')[2])

				sys.exit(1)

			if start > end:
				print "ERROR: Start postion of the target is greater than end position"
				print line
				sys.exit(1)

def fileCheck2(file):
	print "Parsing the target bed file for format check ..."
	
	with open(file) as f:
		for line in f:
			try:
				start = int(line.split('\t')[1])
				end = int(line.split('\t')[2])
			except:
				print "ERROR: Chromosome postions in target bed file are not integers"
				print line
				sys.exit(1)

			if start > end:
				print "ERROR: Start postion of the target is greater than end position"
				print line
				sys.exit(1)				

			try:
				mut = str(line.split('\t')[3].split(';')[0])
			except:
				print "ERROR: Mutation type in target bed file is not snv (Single Nucleotide Variant), ins (Insertion), del (Deletion) or dup (Duplication)"
				print line
				sys.exit(1)

			mutType=[]
			for item in ['snv', 'ins', 'del', 'dup', 'mix']:
				mutType.extend(list(caseIgnore(item)))				
			if not mut in mutType:
				print "ERROR: Mutation type in target bed file is not snv (Single Nucleotide Variant), ins (Insertion), del (Deletion) or dup (Duplication)"
				print line
				sys.exit(1)					

			try:
				maf = float(line.split('\t')[3].split(';')[1])
			except:
				print "ERROR: Mutation fraction in target bed file is not a number between 0-1"
				print line
				sys.exit(1)

			if not 0 <= maf <= 1:
				print "ERROR: Mutation fraction in target bed file is not a number between 0-1"
				print line
				sys.exit(1)

			try:
				alt = str(line.split('\t')[3].split(';')[2])	
			except:
				print "ERROR: Alternate sequence field in target bed file doesn't exist"
				print line
				sys.exit(1)	
			if (not validSequence(alt)):
				print "ERROR: Alternate sequence specified in target bed file contains characters other than A/C/G/T."
				print line
				sys.exit(1)	

			try:
				length = line.split('\t')[3].split(';')[3]
			except:
				print "ERROR: Indel length field in target bed file doesn't exist"
				print line
				sys.exit(1)	

			if type(length) is int:
				if not int(length) > 0:
					print "ERROR: Indel length specified in target bed file is not a positive integer; default length of 10 bp will be used"
					print line
			if ((type(length) is int) and (type(alt) is str)):
				if (int(length) > 0) and (len(alt) > 0):
					if int(length) != len(alt):
						print "WARNING: Length of insertion sequence in does not match with indel length specified in target bed file. Indel length will be ignored and the specified insertion sequence will be used"
						print line


def caseIgnore(s):
    return (''.join(t) for t in itertools.product(*zip(s.lower(), s.upper())))

def validSequence(s):
    valid = 'ACTG'
    for letter in s:
        if letter not in valid:
            return False
    return True

def get_mut_location(bedentry):
	entryitems = bedentry.split('\t')
	entrychrom = entryitems[0]
	mutposition = int(entryitems[1]) + (int(entryitems[2]) - int(entryitems[1]))/2
	return entrychrom + '\t' + str(mutposition)

def insert(sequence, position, insertseq):
	try:
		index_seq=list(sequence)
		sequencelength = len(sequence)
		index_seq.insert(position, insertseq)
		newsequence = ''.join(index_seq)
		return newsequence[0:sequencelength]
	except:
		print "Exception."
		print sequence
		print position
		return sequence

def insertAMP(query, sequence, position, insertseq, genome, mutation_chrom, mutation_start, amplicon, readlength):
	try:
		index_seq=list(sequence)
		sequencelength = len(sequence)
		index_seq.insert(position, insertseq)
		newsequence = ''.join(index_seq)
		amplicon_pos=checkPosBed(mutation_chrom,mutation_start,amplicon)

		if (query.is_reverse == False):
			extrasequence_start = mutation_start + sequencelength - position - query.query_alignment_start + 1
			extrasequence_end = amplicon_pos[2]
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = newsequence + extrasequence + ("N"*readlength)
			print "Forward Sequence: %s" %newsequence[0:readlength]
			return newsequence[0:readlength]

		else:
			extrasequence_start = amplicon_pos[1]
			extrasequence_end = mutation_start - position - query.query_alignment_start
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = ("N"*readlength) + extrasequence + newsequence
			return newsequence[(len(newsequence)-readlength):]				

	except:
		print "Exception."
		print sequence
		print position		
		return sequence


def revcomp(str):
	out = []
	tmp = list(str)
	tmp.reverse()
	for c in tmp:
		if(c.upper() == "A"):
			out.append("T")
		elif(c.upper() == "T"):
			out.append("A")
		elif(c.upper() == "C"):
			out.append("G")
		elif(c.upper() == "G"):
			out.append("C")
		elif(c.upper() == "N"):
			out.append("N")
		else:
			out.append("!")
	outputstring="".join(out)
	return outputstring

	
def reversed_string(a_string):
    return a_string[::-1]


def quality_score_encoder(qscorearray):
	outputqscores = ""
	for item in qscorearray:
		outputqscores = outputqscores + chr(item + 33)
	return outputqscores


def deletionAMP(query,sequence, position, dellength, genome, mutation_chrom, mutation_start,amplicon,readlength):
	try:
		sequencelength = len(sequence)
		delendposition = position + dellength
		amplicon_pos=checkPosBed(mutation_chrom,mutation_start,amplicon)
		if (query.is_reverse == False):
			extrasequence_start = mutation_start + sequencelength - position - query.query_alignment_start + 1
			extrasequence_end = amplicon_pos[2]
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = sequence[0:position] + sequence[delendposition:] + extrasequence + ("N"*readlength)
			return newsequence[0:readlength]
		else:
			extrasequence_start = amplicon_pos[1]
			extrasequence_end = mutation_start - position - query.query_alignment_start
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = ("N"*readlength) + extrasequence + sequence[0:position] + sequence[delendposition:]
			return newsequence[(len(newsequence)-readlength):]				
	except:
		print "Exception."
		print sequence
		print position
		return sequence	


fragments = {}
def deletion(infile, query, genome, mutation_chrom, mutation_start, indellen):
	time_start = datetime.datetime.now().replace(microsecond=0)
	if ((query.reference_id == query.next_reference_id) and (query.is_proper_pair)):
		mate = infile.mate(query)
		if ((query.is_reverse == False) and (query.reference_start <= mate.reference_start)):
			forward = query
			reverse = mate

			# IF QUERY NAME DOESN'T EXIST IN FRAGMENT DICTIONARY -> 	Run frgament function which returns fragment sequence and create new key (read name)/ value (fragment sequence) pair in dictionary 
			if fragments.has_key(query.query_name) == False:
				fragments[query.query_name]	 = fragmentCreator(forward, reverse, genome, mutation_chrom, mutation_start, indellen)
				
			# Retrun first 101 bases from fragment sequence
			return fragments[query.query_name][0:forward.query_alignment_length]

		elif ((query.is_reverse == True) and (query.reference_start >= mate.reference_start)):
			forward = mate
			reverse = query

			#IF QUERY NAME DOESN'T EXIST IN FRAGMENT DICTIONARY -> 	Run frgament function which returns fragment sequence and create new key (read name)/ value (fragment sequence) pair in dictionary 
			if fragments.has_key(query.query_name) == False:
				fragments[query.query_name]	 = fragmentCreator(forward, reverse, genome, mutation_chrom, mutation_start, indellen)				

			#Retrun last 101 bases from fragment sequence
			return fragments[query.query_name][-reverse.query_alignment_length:] 


def fragmentCreator(forward, reverse, genome, mutation_chrom, del_start, dellength):
	del_end = del_start + dellength

	# First generate the fragement using forward and reverse read information
	fragment = forward.query_sequence[forward.query_alignment_start:forward.query_alignment_end]
	if ((reverse.reference_start - (forward.reference_end - 1)) >= 0):
		fragment = fragment + genome[mutation_chrom][(forward.reference_end + 1):(reverse.reference_start + 1)] # get ref seq
		fragment = fragment + reverse.query_sequence[reverse.query_alignment_start:reverse.query_alignment_end]
	else:
		fragment = fragment + reverse.query_sequence[reverse.query_alignment_start:reverse.query_alignment_end][((forward.reference_end - 1) - reverse.reference_start + 1):]

	if del_end <= reverse.reference_end:
		fragment = fragment[:(del_start - forward.reference_start)] + fragment[(del_end - forward.reference_start):] + genome[mutation_chrom][(reverse.reference_end + 1):(reverse.reference_end + dellength + 1)]		
	else:
		fragment = fragment[:(del_start - forward.reference_start)] + genome[mutation_chrom][(del_end + 1):(del_end + dellength + 1)]

	return fragment


def duplication(sequence, position, duplength, genome, mutation_chrom, mutation_start):
	try:
		dupsequence = genome[mutation_chrom][mutation_start+1:mutation_start + duplength+1]
		sequencelength = len(sequence)
		newsequence = sequence[0:position] + dupsequence + sequence[position:]
		return newsequence[0:sequencelength], dupsequence
	except:
		print "Exception."
		print sequence
		print position
		return sequence

def duplicationAMP(query, sequence, position, duplength, genome, mutation_chrom, mutation_start, amplicon, readlength):
	try:
		sequencelength = len(sequence)
		amplicon_pos=checkPosBed(mutation_chrom,mutation_start,amplicon)
		dupsequence_end = amplicon_pos[2]
		dupsequence = genome[mutation_chrom][(mutation_start+1):dupsequence_end]
		dupsequence = dupsequence[:duplength]
		newsequence = sequence[0:position] + dupsequence + sequence[position:]
		
		if (query.is_reverse == False): 
			extrasequence_start = mutation_start + sequencelength - position - query.query_alignment_start + 1
			extrasequence_end = amplicon_pos[2]
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = newsequence + extrasequence + ("N"*readlength)
			print "Forward Sequence: %s" %newsequence[0:readlength]
			return newsequence[0:readlength], dupsequence
		else:
			extrasequence_start = amplicon_pos[1]
			extrasequence_end = mutation_start - position - query.query_alignment_start
			extrasequence = genome[mutation_chrom][extrasequence_start:extrasequence_end]
			newsequence = ("N"*readlength) + extrasequence + newsequence
			return newsequence[(len(newsequence)-readlength):], dupsequence	
	
	except:
		print "Exception."
		print sequence
		print position
		return sequence

def snv(sequence, readposition, genome, chrom, chromposition, alt):
	try:
		genomebase = genome[chrom][chromposition].upper()
		if alt:
			newbase = alt[0]
		elif genomebase == "A":
			newbase = "T"
		elif genomebase == "T":
			newbase = "A"
		elif genomebase == "C":
			newbase = "G"
		elif genomebase == "G":
			newbase = "C"
		else:
			newbase = "N" 
		index_seq=list(sequence)
		index_seq[readposition-1]=newbase
		return ''.join(index_seq), newbase
	
	except:
		print "Exception."
		print sequence
		print readposition
		return sequence


def get_mut_location(bedentry):
	entryitems = bedentry.split('\t')
	entrychrom = entryitems[0]
	mutposition = int(entryitems[1]) + (int(entryitems[2]) - int(entryitems[1]))/2
	return entrychrom + '\t' + str(mutposition)


def readAmpBed(file):
	rev_exl=open(file,'r')
	RX=[]
	for line in rev_exl:
		line=line.rstrip()
		nline=line.split('\t')
		RX.append(nline)
	rev_exl.close()
	return RX


def checkPosBed(this_chrm,this_pos,RX):
	for pos in RX:
		chrm=pos[0]; start=int(pos[1]); stop=int(pos[2])
		if(this_chrm==chrm and this_pos>=start and this_pos<=stop):
			return [chrm,start,stop]
			break
	return False


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-config", dest="config")

	args = parser.parse_args()

	if (args.config == None):
		print
		print "insiM: Genomic In-Silico Mutator Program for Bioinformatics Validation of Clinical Next Generation Sequencing Assays"
		print "Genomic and Molecular Pathology, University of Chicago"
		print "Version 1.0.0 - 8/21/2018"
		print
		print "FAILED: insiM REQUIRES THE FOLLOWING PARAMETERS:"
		print
		print "\t** REQUIRED PARAMETERS **"
		print "\t-config Configuration file with simulation parameters"

		sys.exit(1)

	params = configParser(args.config)
	assay = params["assay"]
	bamfile = params["bam"]
	targetbed = params["target"]
	mutationtype = params["mutation"]
	genomefasta = params["genome"]
	amplicon = params["ampliconbed"]
	readlength = params["read"]
	mutfraction = params["vaf"]
	indellength = params["len"]
	insertsequence = params["seq"]
	outputfastqname = params["out"]


	if (assay == None) or (bamfile == None) or (targetbed == None) or (mutationtype == None) or (genomefasta == None) or (mutationtype not in list(caseIgnore('mix')) and mutfraction == None) or (assay in list(caseIgnore('amplicon')) and (amplicon == None) or (readlength == None)):
		print
		print "insiM: Genomic In-Silico Mutator Program for Bioinformatics Validation of Clinical Next Generation Sequencing Assays"
		print "Genomic and Molecular Pathology, University of Chicago"
		print "Version 1.0.0 - 8/21/2018"
		print
		print "FAILED: insiM CONFIGURATION FILE REQUIRES THE FOLLOWING PARAMETERS:"
		print
		print "\t** REQUIRED PARAMETERS **"
		print "\t-assay ASSAY TYPE: 'amplicon' or 'capture' (Hybrid Cature)"
		print "\t-bam ALIGNED BAM FILE"
		print "\t-target MUTATION TARGET BED FILE"
		print "\t-mutation MUTATION TYPE: 'snv' (Single Nucleotide Variant), 'ins' (Insertion), 'del' (Deletion), 'dup' (Duplication) or 'mix' (Mixed)>"
		print "\t-genome GENOME FASTA FILE"
		print
		print "\t** CONDITIONAL/OPTIONAL PARAMETERS **"
		print "\t-ampliconbed AMPLICON BED FILE. REQUIRED WHEN -assay (ASSAY TYPE) == 'amplicon'"
		print "\t-read READ LENGTH. REQUIRED WHEN -assay (ASSAY TYPE) == 'amplicon'"
		print "\t-vaf MUTATION FRACTION (MUST BE INTEGER, 1/RATE). REQUIRED WHEN -mutation (MUTATION TYPE) != 'mix' (Mixed)"
		print "\t-len SIZE OF INSERT. IF NOT SPECIFIED, RANDOM SEQUENCE OF LENGTH 10 WILL BE USED FOR ins/del/dup"
		print "\t-seq SEQUENCE OF INSERT. FOR INSERTION, THIS SEQUENCE WILL BE INSERTED; IF NOT SPECIFIED, RANDOM SEQUENCE OF LENGTH SPECIFIED BY -len WILL BE INSERTED."
		print "\t-out OUTPUT FASTQ BASENAME - BASENAME_R1_001.fastq."

		sys.exit(1)

	main(assay, bamfile, targetbed, mutfraction, outputfastqname, indellength, insertsequence, mutationtype, genomefasta, amplicon, readlength)
