#! /usr/bin/python
# -*- coding: utf-8 -*-

#  _______ _____            _____  
# |__   __|  __ \     /\   |  __ \ 
#    | |  | |__) |   /  \  | |__) |
#    | |  |  _  /   / /\ \ |  ___/ 
#    | |  | | \ \  / ____ \| |     
#    |_|  |_|  \_\/_/    \_\_|     
#                              
# Trasposable elements annotation pipeline


# TOOLS
#
# La felicità è sovversiva se collettivizzata. Andrea Pazienza.
#
# Fabio Barteri		STRUCTURAL EVOLUTION OF PLANT GENOMES GROUP
#					CENTRE FOR RESEARCH IN AGRICULTURAL GENOMICS
#					CARRER DE LA VALL MORONTA, EDIFICI CRAG, CAMPUS UAB
#					08000 CERDANYOLA DEL VALLÉS, BARCELONA, CATALUNYA
#					fabio.barteri@cragenomica.es

import os, sys
from optparse import OptionParser
from itertools import chain

## 1 CLASSES

## 1.1 Class retroelement()		A class to define the single sequences
##								and their features

class retroelement():			
	def __init__(self):
		self.id = ''
		self.pseudo = ''
		self.seed = ''
		self.clust = ''
		self.length = 0
		self.internal = ''
		self.ltr = []
		self.ltr_lengths = []
		self.meanltrlen = 0
		self.deltaltrlen = 0
		self.ltr_proportion = 0
	
	def fetchinternal(self):		# Fetches the internal region
		scaffold = self.pseudo.split(':')[0]
		threeprime_list = self.ltr[0].split(':')[1].split('-')
		fiveprime_list = self.ltr[1].split(':')[1].split('-')				
		lvals = [int(threeprime_list[0]), int(threeprime_list[1]), int(fiveprime_list[0]), int(fiveprime_list[1])]
		lvals.sort()
		
		a = lvals[0]
		b = lvals[1]		# b and c will represent the internal region
		c = lvals[2]
		d = lvals[3]
		
		self.internal = scaffold+':'+str(b)+'-'+str(c)
		
	def ltrcalc(self):
		self.meanltrlen = (self.ltr_lengths[0]+self.ltr_lengths[1])/2
		self.deltaltrlen = abs(self.ltr_lengths[0]-self.ltr_lengths[1])
		self.ltr_proportion = (self.ltr_lengths[0]+self.ltr_lengths[1])*1.0/self.length


## 1.2 class outfile()			Any single output file is described by
##								this class. This is useful to track all
##								the generated files and to print a datasheet
##								containing all the information about the output

class outfile():
	def __init__(self):
		self.name = ''
		self.pwd = ''
		self.description = ''
		self.highlight = False


## 2 FUNCTIONS

## 2.0 rread():				Simply reads a file importing into a list

def rread(filein):
	a = open(filein,'rU')
	alist = a.read().split('\n')
	a.close()
	return alist

## 2.1 launch_harvest() : 	Launches harvest saving the output in a
##							given folder.

'''
Launch suffixerator:
gt suffixerator -db [genome.fa] -indexname [genome.fa] -tis -suf -lcp -des -ssp -sds -dna

Launch LTR harvest:
LTR harvest launch:
gt ltrharvest -index [genome.fa] -gff3 [OUTPUT.gff3] #StandardOptions
'''

def launch_harvest(genome_file, output_folder):
	
	suffixerator_command = ' '.join(['gt suffixerator -db', genome_file ,'-indexname', genome_file, '-tis -suf -lcp -des -ssp -sds -dna'])
	os.system(suffixerator_command)
	
	raw_harvest_result = output_folder + '/' + genome_file.split('.')[0]+'.rawannotation.gff3'
	raw_fasta_result = output_folder + '/' + genome_file.split('.')[0]+'.rawfasta.fa'
	harvest_command = ' '.join(['gt ltrharvest -index', genome_file, '-gff3', raw_harvest_result, '-out', raw_fasta_result,'>',output_folder+'/harvest.log'])

	os.system(harvest_command)


## 2.2 harvest_dictionary():	Converts the harvest annotation (0 based)
##								Back to the original annotation
'''
LTR Harvest has the terrible habit to provide the results in a new
annotation, with the prefix Seq and 0-based numeration. This has
awful consequences when running in fragmented genomes, where the
correspondance between number of the scaffold and actual order in the
fasta is not linear (e.g. S. commersonii 2015 assembly).
For this, we need to build a specific dictionary to convert the harvest
made chromosome position back to the normal.

The fasta output will be the only one reporting the conversion between
harvest- generated sequences and the original ones. 
'''

def harvest_dictionary(output_folder, rawfasta, genome_file):
	
	hdict = {}
	
	raw_fasta_headers = output_folder + '/' + genome_file.split('.')[0]+'.raw_fasta_headers.fa'
	id_conversion_file = output_folder + '/' + genome_file.split('.')[0]+'.idconversion.tab'
	
	os.system("grep '>' " + rawfasta + " > " + raw_fasta_headers)
	
	headers = rread(raw_fasta_headers)
	# Headers are in the form: >C2858846_17.0 (dbseq-nr 63598) [5156,11634]
	for header in headers:
		try:
			scaffold = header.split(' ')[0].replace('>','')
			harvest_id = header.split('(')[1].split(')')[0].replace('dbseq-nr ','seq')
			hdict[harvest_id] = scaffold
		except:
			print header
			pass
	
	return hdict
	
## 2.3 filter_harvest_line():		Filters harvest raw result by 

def filter_harvest(hdict, rawfasta, modified_fasta):
	handle = rread(rawfasta)
	outhand = open(modified_fasta, 'w')
	for l in handle:
		try:
			c = l.split('\t')
			harvest_sequence = c[0]
			real_sequence = hdict[harvest_sequence]
			print >> outhand, l.replace(harvest_sequence, real_sequence) 
		except:
			print >> outhand, l
			pass
	outhand.close()

## 2.4 bedline():					From ID and pseudoid generates a bedlilne

def bedline(id, pseudoid, score = '.', strand = '.'):
	scaffold = pseudoid.split(':')[0]
	start = pseudoid.split(':')[1].split('-')[0]
	stop = pseudoid.split(':')[1].split('-')[1]
	return '\t'.join([scaffold, start, stop, id, score, strand])


## 2.5 getfasta():					It is used to fetch the fasta from a gff/bed
#									a simple bedtools getfasta workaround

def getfasta(bed_file, genome_file, fasta_outfile, name = 'yes', strand = 'no', more_options = ''):
	bt_head = 'bedtools getfasta'
	bed_input = ' -bed ' + bed_file
	fasta_input = ' -fi ' + genome_file
	fasta_output = ' -fo ' + fasta_outfile
	
	options = ' '
	if name == 'yes':
		options = options + '-name '
	if strand =='yes':
		options = options + '-s'
	options = options + more_options
	
	command = bt_head + options + fasta_input + bed_input + fasta_output
	os.system(command)

## 2.6 launch_sixframe():			This will launch sixframe translation tool
#									used for hmmer - based domain scanning

def launch_sixframe(fasta_in, sixfasta_out):
	# transeq -sequence (dna.fa) -outseq (six_frame.fa) -frame 6
	command = ' '.join(['transeq -sequence', fasta_in, '-outseq', sixfasta_out, '-frame 6'])
	os.system(command)

## 2.7 lanuch_hmmscan():			The function to launch hmmscan
	
def launch_hmmscan(domain_table, sequence_table, profiles, sixframe,hlog):
	# hmmscan ­­--domtblout SixFrame_per_domain.tsv ­­--tblout SixFrame_per_sequence.tsv ­­ --noali ­ GyDb_profiles.hmm SixFrame_prot.fa
	command = ' '.join(['hmmscan --domtblout', domain_table, '--tblout', sequence_table, '--noali', profiles, sixframe,'>',hlog])
	os.system(command)

## 2.8 lanuch_trf():			The function to launch trf internal repetitions finder

def launch_trf(trf_folder, dat_file, fasta_infile, logfile):
	#trf yoursequence.txt 2 7 7 80 10 50 500 -f -d 
	params = '2 7 7 80 10 50 500 -f -d'
	trf_command = ' '.join(['trf', fasta_infile, params, '>', logfile])
	html_move_command = ' '.join(['mv', '*html',trf_folder])
	dat_cat_command = ' '.join(['cp', '*dat',dat_file])
	os.system(trf_command)
	os.system(html_move_command)
	os.system(dat_cat_command)
	os.system('rm -r *dat')

## 2.9 Join Ranges (from the internet, I need to study this, but it works!)

flatten = chain.from_iterable

LEFT, RIGHT = 1, -1

def join_ranges(data, offset=0):
    data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT))
            for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield x, value - offset

## 2.10 sum_intervals			Given a list of interval tuples [(start1, stop1), (start2, stop2),...]
##								it returns the sum of the intervals. This is to calculate the trf_entry.repsum
##								from the trf_entry.merged_repintervals.

def sum_intervals(list_of_tuples):
	s = 0
	for x in list_of_tuples:
		su = abs(x[1]-x[0])
		s = s + su
	
	return s


## 2.11 process_trf():			The function to process trf output and generate the conversion table
def process_trf(dat_file, outuput_table, length_dictionary):
	
	# Step 0: define the trf_entry class
	class trf_entry():
		
		def __init__(self):
			self.id = ''					# The element ID
			self.element_length = 0			# The element length
			self.repintervals = []			# The list of repetitive intervals [(start1, stop1), (start2, stop2),...] 
			self.merged_repintervals = []	# The list of merged repetitive intervals (they overlap and need to be merged)
			self.repsum = 0.0				# The sum of the lengths of repetitive interval. It gives the total nts that are tandem repetitions.
			self.ratio = 0.0				# The ratio between repsum/element_length
			self.tag = 'NO_INTERNAL_REPS'	# The final tag (NO_INTERNAL_REPS for ratio < .5, AUTOREPETITIVE.flag[ratio] for ratio > .5) 
	
		def loadinfo(self):
			self.merged_repintervals = list(join_ranges(self.repintervals))
			self.repsum = sum_intervals(self.merged_repintervals) * 1.0 # has to be float
			self.ratio = self.repsum/self.element_length
			
			if self.ratio > 0.5:
				self.tag = 'AUTOREPETITIVE-flag'+ str(self.ratio)[1:4]
			
	
	#####	Start calculation	#####
			
	# Step 1 : import the dat_file and parse it for sequence
	infile = open(dat_file, 'rU')
	inhand = infile.read()
	infile.close()
	
	parsed_for_sequence = inhand.split('Sequence: ')
	
	# Step 2: import the autorepetitive intervals for each sequence
	
	trfentries = []
	
	for entry in parsed_for_sequence:
		try:
			splitlist = entry.split('\n')
			z = trf_entry()
			z.id = splitlist[0]
			z.element_length = length_dictionary[z.id]

			for l in splitlist:
				if 'Parameters' not in l:
					try:
						c = l.split(' ')
						interval = [int(c[0]), int(c[1])]
						z.repintervals.append(interval)
					except:
						pass
			z.repintervals.sort()
			trfentries.append(z)
		except:
			pass
	
	# Step 3: process the auto-repetitiveness and print the output_table
	
	outhand = open(outuput_table,'w')
	
	for z in trfentries:
		z.loadinfo()
		print>>outhand, '\t'.join([z.id, z.tag])
	
	outhand.close()

## 2.12 hasher():			From a tabbed file, it returns a dictionary
##							where the key is the first column and the value
##							is the second


def hasher(filein, mode = 'string'):
	outdict = {}
	infile = open(filein, 'rU')
	inlist = infile.read().split('\n')
	infile.close()
	
	for l in inlist:
		try:
			c = l.split('\t')
			key = c[0]
			value = c[1]
			if mode == 'numeric':
				value = int(value)
			
			outdict[key] = value
		except:
			pass
	
	return outdict


## 2.13 launch_makeblastdb():	Creates a blast database, requires NCBI blast

def launch_makeblastdb(genome_file):
	#makeblastdb -in peach_reduced_100k.fa -input_type fasta -dbtype nucl
	command = ' '.join(['makeblastdb -in', genome_file, '-input_type fasta -dbtype nucl'])
	os.system(command)
			

## 2.14 launch_blast():			Launch blast (outfmt = 6, tabbed file as required by SiLiX)

def launch_blast(infile, dbname, outfile):
	command = ' '.join(['./blastn', '-query', infile, '-db', dbname, '-outfmt 6', '-out',outfile])
	os.system(command)

## 2.15 launch_silix(): 		Launch SiLiX

def launch_silix(fastafile, blastfile, prefix, outfile):
	# silix [OPTIONS] <FASTAFILE> <BLASTFILE>
	options = '-i 0.8 -r 0.8 -f'
	command = ' '.join(['silix', options, prefix, fastafile, blastfile, '>', outfile])
	os.system(command)

## 2.16 process_silix():		Process the silix result

def process_silix(fnodes_file, id2clust):
	clusters = []
	clust2num = {}
	clust2tag = {}
	tag = 'SINGLE'
	
	inhand = open(fnodes_file, 'rU')
	inlist = inhand.read().split('\n')
	inhand.close()
	
	outhand = open(id2clust,'w')
	# Fetch the clusters
	
	for l in inlist:
		try:
			c = l.split('\t')
			clust_id = c[0]
			ltr_id = c[1]
			
			if clust_id not in clusters:
				clusters.append(clust_id)
		except:
			pass
	clusters.sort()

	# Count the members and sort between clustered and singletons
	for cl in clusters:
		count = 0
		for l in inlist:
			try:
				c = l.split('\t')
				clust_id = c[0]
				
				if cl == clust_id:
					count = count + 1
				
			except:
				pass
		
		clust2num[cl] = str(count)
		
		if count > 1:
			tag = 'CLSTR'
		else:
			tag = 'SINGLE'
		
		
		clust2tag[cl] = tag
	
	# And print the table

	for l in inlist:
		try:
			c = l.split('\t')
			clust_id = c[0]
			ltr_id = c[1]
			
			print >> outhand, '\t'.join([ltr_id, clust_id, clust2num[clust_id], clust2tag[clust_id]])
		except:
			pass

			
			
			

