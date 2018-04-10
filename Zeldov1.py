#Python-Zeldovich v1


from Bio import SeqIO
import itertools
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob

#use SeqIO module from Biopthyon
IVYWREL = 0
#create and set the counter varable for IVYWREL amino acids
NONE = 0
#create and set the counter variable for other amino acids
GAP = 0
COUNTER = 0

file = open ("results2.txt","w")
#this opens the the text file 'results.txt', for writing
file.write (' Sequence counter: O')
files = glob.glob('*.phy')
for x in files:
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
file.close()


#this closes the file as we'd opened it for writing
with open('results2.txt', 'r') as f:
	OGTLINES = itertools.islice(f, 5, None, 6)
	contentOGT = [x.strip() for x in OGTLINES]
	#print contentOGT
	contentOGTfloatlist = map(float, contentOGT)
	totalcontentOGT = sum(contentOGTfloatlist)
	#print totalcontentOGT
	numberofsequences = len(contentOGT)
	#print numberofsequences
	MEANOGT =totalcontentOGT / numberofsequences
	#print "MEAN OGT:  ", MEANOGT
	MINOGT = min(contentOGTfloatlist)
	#print "MIN OGT: ", MINOGT
	MAXOGT = max(contentOGTfloatlist)
	#print "MAX OGT:  ", MAXOGT
	RANGE = MAXOGT - MINOGT
	#print 'RANGE:  ', RANGE
with open ('results2.txt', 'r') as f:
	FRACTION = itertools.islice(f, 4, None, 6)
	contentFRACTION = [x.strip() for x in FRACTION]
	#print 'FRACTION OF IVYWREL: ', contentFRACTION

#plt.plot([contentFRACTION], [(937 * contentFRACTION) - 355])
#temperaturelist = sum[(937 * contentFRACTION) - 335]
plt.plot(contentOGTfloatlist, contentFRACTION)
plt.ylabel ('F of IVYWREL Amino Acids')
plt.xlabel ('Temperature in celsius')
plt.show()



print numberofsequences
print "MEAN OGT: ", MEANOGT
print "MIN OGT: ", MINOGT
print "MAX OGT: ", MAXOGT
print "RANGE: ", RANGE