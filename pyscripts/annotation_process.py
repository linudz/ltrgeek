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


# CORE SCRIPT: THE BODY OF THE PIPELINE
#


# FEB 13, 2017 - MODIFIED AND LIMITED TO LTR HARVEST RUN AND PROCESSING


# La felicità è sovversiva se collettivizzata. Andrea Pazienza.
#
# Fabio Barteri		STRUCTURAL EVOLUTION OF PLANT GENOMES GROUP
#					CENTRE FOR RESEARCH IN AGRICULTURAL GENOMICS
#					CARRER DE LA VALL MORONTA, EDIFICI CRAG, CAMPUS UAB
#					08000 CERDANYOLA DEL VALLÉS, BARCELONA, CATALUNYA

#					fabio.barteri@cragenomica.es

'''
Introduction to TRAP
----------------------------------------------------------------
This pipeline is thought to characterise the LTR transposable elements
from raw genomic data. A genome is scanned with LTR harvest for a
structural identification of putative LTR retroelements. Then,
characterisation goes along with filtering of sequences.

- Sequences are filtered for their proportion of Ns (more than 10%)
- Protein domains and strand orientation are identified with HMMscan
- Insertion time is determined with KIMURA algorithm
- 

'''

from trap_tools import *

parser = OptionParser()
parser.add_option("-g", "--genome", dest="genome", help="The input genome. Warning: the pipeline will take all the sequences into account. If you want to eliminate some (e.g. unassembled scaffolds), please remove them first.")
parser.add_option("-j", "--jobfolder", dest="jobfolder", default="job", help= 'Assign a name to your run. This will generate a folder containing results, logs and intermediate files.')
parser.add_option("-s", "--specie", dest="specie", default="Ppersica", help= 'The name of the specie. This will name the output files')
parser.add_option("-k", "--mutationrate", dest="mutation_rate", default="1.8E-8", help= 'The mutation rate to submit to kimura algorithm for insertion time determination')

(options, args) = parser.parse_args()

# STEP 0: Create the job folder and subfolder

job_folder = options.jobfolder
harvest_folder = job_folder + '/harvest_data'
hmmscan_folder = job_folder + '/hmmscan'
age_folder = job_folder + '/date_divergence'
age_folder_bed = age_folder + '/bedfiles'
age_folder_fasta = age_folder + '/fastafiles'
N_folder = job_folder + '/Ncount'
silix_folder = job_folder + '/SiLiX'
trf_folder = job_folder + '/InternalRepetitions'
highlights = job_folder + '/highlights'

try:
	# To plant the tree
	os.mkdir(job_folder)
	os.mkdir(harvest_folder)
	os.mkdir(hmmscan_folder)
	os.mkdir(age_folder)
	os.mkdir(age_folder_bed)
	os.mkdir(age_folder_fasta)
	os.mkdir(N_folder)
	os.mkdir(silix_folder)
	os.mkdir(trf_folder)
	os.mkdir(highlights)
except:
	# Remove the previous folder and subdirs
	os.system('rm -r ' + job_folder)
	
	# Re-plant the tree
	os.mkdir(job_folder)
	os.mkdir(harvest_folder)
	os.mkdir(hmmscan_folder)
	os.mkdir(age_folder)
	os.mkdir(age_folder_bed)
	os.mkdir(age_folder_fasta)
	os.mkdir(N_folder)
	os.mkdir(silix_folder)
	os.mkdir(trf_folder)
	os.mkdir(highlights)


# Create symbolic link to the genome file
genome = options.genome
#genome_link = 'ln -s ../' + genome
#os.system(genome_link)

filelist = []

##	STEP 1: LTR HARVEST RUN
'''
LTR harvest is first ran starting from the genome. The result is modified
to match with the genome naming of the scaffold (see trap_tools.py 2.2).
'''
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################


##		1.1 Launch harvest

print '[RUNNING LTR HARVEST]'
launch_harvest(genome, harvest_folder) # trap_tools.py - 2.1

##		1.2 Correct the harvest ids back to normal scaffolds/chromosomes

print '[PROCESING RAW LTR HARVEST RESULT] Id correction'

###!! OUTFILE: LTR harvest raw fasta 
raw_fasta_result = outfile()
raw_fasta_result.name = 'LTR harvest raw fasta'
raw_fasta_result.pwd = harvest_folder + '/' + genome.split('.')[0]+'.rawfasta.fa'
raw_fasta_result.description = 'The raw ltr harvest fasta output. This is needed for id conversion and hmmscan exec.'
filelist.append(raw_fasta_result)
#############################################!! OUTFILE !!###

###!! OUTFILE: LTR harvest raw gff
raw_harvest_gff = outfile()
raw_harvest_gff.name = 'LTR harvest raw gff'
raw_harvest_gff.pwd = harvest_folder + '/' + genome.split('.')[0]+'.rawannotation.gff3'
raw_harvest_gff.description = 'The raw output of LTR harvest in gff format. Is the base of any further analysis.'
filelist.append(raw_harvest_gff)
#############################################!! OUTFILE !!###

###!! OUTFILE: Corrected GFF
modified_gff = outfile()
modified_gff.name = 'LTR harvest modified gff'
modified_gff.pwd = harvest_folder + '/' + genome.split('.')[0]+'.corrected_annotation.gff3'
modified_gff.description = 'The LTR harvest gff output with the scaffold/chr name corrected'
filelist.append(modified_gff)
#############################################!! OUTFILE !!###

harvestid2scaffold = harvest_dictionary(harvest_folder, raw_fasta_result.pwd, genome)
filter_harvest(harvestid2scaffold, raw_harvest_gff.pwd, modified_gff.pwd)


##	STEP 2: HARVEST PROCESSING (former hs_core)
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################

##		2.1 Reading the harvest corrected gff

harv = open(modified_gff.pwd,'rU')
inextlist = harv.read().split('\n###\n')
harv.close()

##		2.2 Fetching the single element information
######################################################

print '[PROCESING RAW LTR HARVEST RESULT] Fetching single element information'

el_list = [] # the list of retroelement() objects
ltr_ends_filelist = []		# list of the outfile() objects relative to ltr_ends
							# this is used in date_divergence calculation

for el in inextlist:
	try:
		z = retroelement()
		lines = el.split('\n')
		
		## 2.2.1 : Fetch the element ID, pseudoid, length
		for l in lines:
			c = l.split('\t')
			if c[2] == 'LTR_retrotransposon':
				scaffold = c[0]
				start = int(c[3])
				stop = int(c[4])
				z.length = abs(stop-start)
				z.pseudo = scaffold+':'+str(start)+'-'+str(stop)
				z.id = c[8].split(';')[0].replace('ID=','')
				
		
		## 2.2.2: Fetch the LTR ends in form of pseudoids
		for l in lines:
			c = l.split('\t')
			if c[2] == 'long_terminal_repeat':
				scaffold = c[0]
				start = c[3]
				stop = c[4]
				length = abs(int(stop) - int(start))
				z.ltr_lengths.append(length)
				ltr_pseudo = scaffold+':'+start+'-'+stop
				z.ltr.append(ltr_pseudo)		

		z.ltr.sort()
		z.fetchinternal()
		z.ltrcalc()
		el_list.append(z)
	except:
		pass

## 		2.3 ID <-> PseudoID conversion table generation
######################################################

###!! OUTFILE: ID <-> PseudoID conversion table
id2pseudofile = outfile()			
id2pseudofile.name = 'id2pseudo-convtable.tab'
id2pseudofile.pwd = harvest_folder + '/' + id2pseudofile.name
id2pseudofile.description = 'File to convert ID to PSEUDOID'
id2pseudofile.highlight = True
filelist.append(id2pseudofile)
#############################################!! OUTFILE !!###

print '[PROCESING RAW LTR HARVEST RESULT] ID2position conversion table'
out_hand = open(id2pseudofile.pwd,'w')

for z in el_list:
	print >> out_hand, z.id + '\t' + z.pseudo

out_hand.close()

## 		2.4 ID <-> Element length conversion table generation
######################################################

###!! OUTFILE: ID <-> PseudoID conversion table
id2lenfile = outfile()			
id2lenfile.name = 'id2length.tab'
id2lenfile.pwd = harvest_folder + '/' + id2lenfile.name
id2lenfile.description = 'File to compare id and length'
id2lenfile.highlight = True
filelist.append(id2lenfile)
#############################################!! OUTFILE !!###

print '[PROCESING RAW LTR HARVEST RESULT] ID <-> Element length conversion table generation'


out_hand = open(id2lenfile.pwd,'w')

for z in el_list:
	print>> out_hand, z.id + '\t' + str(z.length)

out_hand.close()


## 		2.5 ID <-> pseudo <-> internal conversion table generation
######################################################

print '[PROCESING RAW LTR HARVEST RESULT] ID <-> pseudo <-> internal conversion table generation'


###!! OUTFILE: ID <-> pseudo <-> internal conversion table
idpsint = outfile()			# File information for id2length.tab
idpsint.name = 'id2pseudo2internal.tab'
idpsint.pwd = harvest_folder + '/' + idpsint.name
idpsint.description = 'File to compare id, pseudoid and internal region pseudoid.\nRequired by convdict.py to read the HMMer results'
idpsint.highlight = True
filelist.append(idpsint)
#############################################!! OUTFILE !!###

out_hand = open(idpsint.pwd,'w')
for z in el_list:
	print>> out_hand, z.id + '\t' + z.pseudo + '\t' + z.internal
out_hand.close()


## 		2.6 Total and complete elements bed file generation
######################################################
print '[PROCESING RAW LTR HARVEST RESULT] Total and complete elements bed and fasta files generation'

###!! OUTFILE: Complete elements bedfile
complete_elements_bed = outfile()
complete_elements_bed.name = options.specie+'-comp_elements.bed'
complete_elements_bed.pwd = harvest_folder + '/' + complete_elements_bed.name
complete_elements_bed.description = 'COMPLETE ELEMENTS Bed file'
filelist.append(complete_elements_bed)
#############################################!! OUTFILE !!###

###!! OUTFILE: Complete elements fastafile
complete_elements_fa = outfile()
complete_elements_fa.name = options.specie+'-comp_elements.fa'
complete_elements_fa.pwd = harvest_folder + '/' + complete_elements_fa.name
complete_elements_fa.description = 'COMPLETE ELEMENTS Fasta file'
filelist.append(complete_elements_fa)
#############################################!! OUTFILE !!###

##		2.6.1 Total and complete elements bed file
out_hand = open(complete_elements_bed.pwd,'w')

for z in el_list:
	l = bedline(z.id, z.pseudo)
	print>> out_hand, l

out_hand.close()

##		2.6.2 Total and complete elements fasta file
getfasta(complete_elements_bed.pwd, options.genome, complete_elements_fa.pwd)


## 		2.7 Internal region bed file
######################################################
print '[PROCESING RAW LTR HARVEST RESULT] Internal region bed and fasta files generation'

###!! OUTFILE: Internal region BEDFILE
intreg = outfile()
intreg.name = options.specie+'-internal_region.bed'
intreg.pwd = harvest_folder + '/' + intreg.name
intreg.description = 'Internal Region Bed file'
filelist.append(intreg)
#############################################!! OUTFILE !!###

###!! OUTFILE: Internal region FASTAFILE (HMMer folder)
intregfa = outfile()
intregfa.name = options.specie+'-internal_region.fa'
intregfa.pwd = hmmscan_folder + '/' + intregfa.name
intregfa.description = 'Internal Region Bed file'
filelist.append(intregfa)
#############################################!! OUTFILE !!###

out_hand = open(intreg.pwd,'w')

for z in el_list:
	l = bedline(z.id, z.internal)
	print>> out_hand, l

out_hand.close()

##			2.7.1 Internal region fasta file (copy in hmmer folder)

getfasta(intreg.pwd, options.genome, intregfa.pwd)

## 		2.8 Single LTR end file generation in the folder
######################################################

print '[PROCESING RAW LTR HARVEST RESULT] Single LTR end file generation for age determination.'

for z in el_list:
	
	###!! OUTFILE(s): LTR ends bedfiles
	bf = outfile()
	bf.name = z.id+'.bed'
	bf.pwd = age_folder_bed + '/' + bf.name
	bf.description = z.id + ' ltr ends bedfile'
	#############################################!! OUTFILE !!###

	###!! OUTFILE(s): LTR ends fastafiles
	fa = outfile()
	fa.name = z.id+'.fa'
	fa.pwd = age_folder_fasta + '/' + fa.name
	fa.description = z.id + ' ltr ends bedfile'
	#############################################!! OUTFILE !!###


	out_hand = open(bf.pwd,'w')
	print>>out_hand, bedline(z.id+'-5ltr', z.ltr[0])
	print>>out_hand, bedline(z.id+'-3ltr', z.ltr[1])
	out_hand.close()
	getfasta(bf.pwd, options.genome, fa.pwd)
	ltr_ends_filelist.append(fa)

