#Python-Zeldovich v1


from Bio import SeqIO
import itertools
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import pandas as pd
import csv
import re

#id = str(386632499)

#species_to_id = {}
#id_to_species = {}

#f = open('nameid.csv')
#for line in f:
#  fields = re.split(",", line.rstrip())
#  species_to_id[fields[0]] = fields[1]
#  id_to_species[fields[1]] = fields[0]

#print species_to_id['Bacillus_anthracis_str_Ames']
#print id_to_species['386632499']

#use SeqIO module from Biopthyon
IVYWREL = 0
#create and set the counter varable for IVYWREL amino acids
NONE = 0
#create and set the counter variable for other amino acids
GAP = 0
COUNTER = 0
ALIGNMENTCOUNTER = 0
TOTALFORALIGNMEAN = 0
SEQCOUNTERPERALIGNMENT = 0
ALIGNMEDTLIST = 0

file = open ("results2.txt","r+")

#this opens the the text file 'results.txt', for writing
files = glob.glob('*.phy')
for x in files: 
	SEQCOUNTERPERALIGNMENT = 0
	ALIGNMENTCOUNTER = ALIGNMENTCOUNTER +1
	ALIGNOGTLIST = []
	file.write ('\n ------------------------NEW ALIGNMENT %s--------------------' %ALIGNMENTCOUNTER)
	FILE = x
	for record in SeqIO.parse (FILE, "phylip"):
		#takes filename and format name and returns a seqrecord  iterator
		#print (record.id)
		sequence=record.seq
		#assign the record to the sequence variable
		for character in sequence:
			#for loop, for each chaaracter in the sequence generated above
				if character in ['I', 'V', 'Y', 'W', 'R', 'E', 'L']:
					#if the characters in the sequence are either I,V,Y,W,R,E or L 
					IVYWREL = IVYWREL + 1
					#increase the counter by one
					# print IVYWREL
				elif character in '-':
					GAP = GAP + 1
				else:
					NONE = NONE + 1
					#otherwise increase the other counter
					# print NONE
		file.write('\n IVYWREL: %s' %IVYWREL)
		file.write(' \n NONE: %s' %NONE)
		#print both counters to the results file

		TOTAL = NONE + IVYWREL
		file.write ('\n TOTAL: %s' %TOTAL)
		FIVYWREL = 1. * IVYWREL / TOTAL
		file.write ('\n %s' %FIVYWREL)
		#this finds the fraction of ivywrel amino acids compared to the other amino acids (excluding gaps)

		OGT = (937 * FIVYWREL) - 335
		TOTALFORALIGNMEAN = TOTALFORALIGNMEAN + OGT
		ALIGNOGTLIST.append (OGT)
		#print ALIGNOGTLIST
		file.write ('\n %s' %OGT) 
		#this uses zeldovich equation to calculate OGT and then prints the output to a file

		IVYWREL = 0
		NONE = 0
		GAP = 0
		TOTAL = 0
		FIVYWREL = 0
		COUNTER = COUNTER + 1
		#this resuts most of the variables used, except COUNTER which is increased by 1 for each loop, to keep track of where we are in the alignment file
		file.write ('\n Sequence counter:  %s' %COUNTER)
		#this prints the sequence counter to the results file
		SEQCOUNTERPERALIGNMENT = COUNTER
		idx = record.id
		file.write ('\n %s' %idx)
		id = str(idx)

		species_to_id = {}
		id_to_species = {}

		f = open('nameid.csv')
		for line in f:
		  fields = re.split(",", line.rstrip())
		  species_to_id[fields[0]] = fields[1]
		  id_to_species[fields[1]] = fields[0]

		#print species_to_id['Bacillus_anthracis_str_Ames']
		file.write (id_to_species[idx])
		
		#file.write ('\n %s' %idx)

	file.write ("\n ----------ALIGNMENT STATISTICS FROM ALIGNMENT %s-----------" %ALIGNMENTCOUNTER)
	ALIGNMEAN = TOTALFORALIGNMEAN / SEQCOUNTERPERALIGNMENT
	#meancalc
	
	#ALIGNMEDLIST = itertools.islice(x, 5, None, 6)	
	#print ALIGNMEDLIST
	#contentalignmed = [g.strip() for g in ALIGNMEDLIST]
	#print contentalignmed
	ALIGNOGTLIST1 = map(float, ALIGNOGTLIST)
	#print ALIGNOGTLIST1

	data = np.array(ALIGNOGTLIST1)
	#print data
	sorted(data)
	ALIGNMED = np.median(data)
	#mediancalc
	min4range = min(data) 
	max4range = max(data)

	ALIGNRANGE = max4range - min4range
	#rangecalc
	#file = open ("alignmentstatfile.txt", "r+")
	#file.write ("\n SEQCOUNTPERALIGN: %s" %SEQCOUNTERPERALIGNMENT)
	#THIS TESTS TO SEE IF THE MEAN IS BEING CALCULATED FROM THE RIGHT TOTAL
	file.write ("\n ALIGNMENT MEAN OGT: %s" %ALIGNMEAN)
	file.write ("\n ALIGNMENT MEDIAN: %s" %ALIGNMED)
	file.write ("\n ALIGNMENT RANGE: %s" %ALIGNRANGE)
	#file.close()
	#file = open ("results2.txt","r+")


file.close()


#this closes the file as we'd opened it for writing

#THE FOLLOWING CODE ANALYSES ALL THE DATA, EVERY SINGLE ALIGNMENT
with open('results2.txt', 'r+') as f:
	OGTLINES = itertools.islice(f, 5, None, 6)
	contentOGT = [x.strip() for x in OGTLINES]
	contentOGTfloatlist = map(float, contentOGT)
	totalcontentOGT = sum(contentOGTfloatlist)
	numberofsequences = len(contentOGT)
	#Looks for the OGT data in the total results file, then assigns the results to a list, the list is then used to create a total for the mean, the number of objects in the content OGT list is used to find out how many sequences there are

	MEANOGT =totalcontentOGT / numberofsequences
	MINOGT = min(contentOGTfloatlist)
	MAXOGT = max(contentOGTfloatlist)
	RANGE = MAXOGT - MINOGT
	#calculates the stat data for all the sequences in all the alignments

with open ('results2.txt', 'r') as f:
	FRACTION = itertools.islice(f, 5, None, 6)
	contentFRACTION = [x.strip() for x in FRACTION]
	#Looks for the fraction of IVYWREL amino acids in the whole results file, and assigns that number to an object

	plt.plot(contentOGTfloatlist, contentFRACTION)
	plt.ylabel ('F of IVYWREL Amino Acids')
	plt.xlabel ('Temperature in celsius')
	plt.show()
	#plots a graph based on all the data

	print 'Number of sequences: ', numberofsequences
	print "MEAN OGT: ", MEANOGT
	print "MIN OGT: ", MINOGT
	print "MAX OGT: ", MAXOGT
	print "RANGE: ", RANGE
	#Prints the stat data to the terminal

#END OF TOTAL ALIGNMENT ANALYSIS
