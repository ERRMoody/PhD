from Bio import SeqIO
from Bio import AlignIO
from Bio import AlignIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

counterofstrains = 0
counterofstrains2 = 0 

def aminoacidchecker(inputfile):
	#this defines a new function, named aminoacidchecker, which takes one input argumnent 'inputfile'



	AAstring = ("GALMFWKQESPVICYHRNDT")
	#This is a string of the 20 aminos acids


	#this is an example string of amino acids to test the code

	#for record in inputfile:
		#print(record)

	#this prints the character for each present in the string of example AAs

	aa_count = {}
	totalcount = 0
	#this creates an empty dictionary

	for aminoacid in AAstring:
		#for each character in the string of all amino acids

		aa_count[aminoacid] = 0

		#the amino acid count variable 20 amino acids is 0

	for record in inputfile:
		if record in aa_count:
			aa_count[record] += 1
			totalcount = totalcount + 1

		else:
			exit
	#print aa_count[aa]
		
	for character in AAstring:
		print ("\n --- " + character + " ---")
		print ("There are " + str(aa_count[character]) + " copies of the amino acid " + str(character) + ", out of a total of " + str(totalcount))
		foo = aa_count[character]
		#print (foo)
		bar = totalcount
		#print (bar)
		percentageofAA = float(foo)/float(bar)
		#percentageofAA = foo/bar
		percentageofAA = percentageofAA 	
		print("Proportion of " + character + " present is " + str(percentageofAA) + "")



alignment = AlignIO.read(open("6.phy"), "phylip")
print("Alignment length %i" % alignment.get_alignment_length())


for record in alignment:
	#print(record.seq + " " + record.id)
	counterofstrains2 = counterofstrains2 + 1
	print ("\n\n\n-------------------------------- New strain " + str(counterofstrains2) + "---------------------------------\n")

	aminoacidchecker(record)

	beebop = record.id

	print ("\n -------------------------------- Strain ID  " + beebop + " --------------------\n")


print ("\n\n\n ########################################## \n\n\n")

alignment = AlignIO.read(open("3927.phy"), "phylip")
print("Alignment length %i" % alignment.get_alignment_length())


for record in alignment:
	#print(record.seq + " " + record.id)
	counterofstrains = counterofstrains + 1
	print ("\n\n\n-------------------------------- New strain " + str(counterofstrains) + "---------------------------------\n")

	aminoacidchecker(record)

	beebop = record.id

	print ("\n -------------------------------- Strain ID  " + beebop + " --------------------\n")



#Comment below is code to convert a pandasdataframe to a csv file
#DataFrame.to_csv(path_or_buf=None, sep=', ', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', encoding=None, compression=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, tupleize_cols=None, date_format=None, doublequote=True, escapechar=None, decimal='.')
