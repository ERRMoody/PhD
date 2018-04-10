import Bio
from Bio import SeqIO
from Bio import AlignIO
import itertools
import numpy as np
import glob
import pandas as pd
import csv
import re
from collections import defaultdict
#Packages being imported

dipsy = []
species_to_id = {}
id_to_species = {}

species_ogts = defaultdict(list)

f = open('nameid.csv')
for line in f:
	fields = re.split(",", line.rstrip())
	species_to_id[fields[1]] = fields[2]
	id_to_species[fields[2]] = fields[1]

#print species_to_id['Bacillus_anthracis_str_Ames']

#a reversible species dictionary is set up to link localIDs in the text to the strains used in the alignments

ALIGNMENTCOUNTER = 0
COUNTER = 0
GAP = 0
SEQCOUNTERPERALIGNMENT = 0
AMINOACIDSPERSTRAIN = 0
NONE = 0

G = 0
A = 0
L = 0
M = 0
F = 0
W = 0
K = 0
Q = 0
E = 0
S = 0
P = 0
V = 0
I = 0
C = 0
Y = 0
H = 0
R = 0
N = 0
D = 0
T = 0

#This is a list of all amino acids, these will evenetually become percentages, but are being created here for use later on


file = open ("fixed4.csv", "w+")
file.write ('Strain,Local ID,Alignment file,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y')

#this opens a results file for writing
files = glob.glob('*.phy')
for x in files:
	SEQCOUNTERPERALIGNMENT = 0
	ALIGNMENTCOUNTER = ALIGNMENTCOUNTER +1
	AMINOACIDSPERSTRAIN = 0
	#print ('\n ALIGNMENT COUNTER %s' %ALIGNMENTCOUNTER)
	FILE = x
	for record in SeqIO.parse (FILE, "phylip"):
		sequence=record.seq
		for character in sequence:
			if character in ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T']:
				SEQCOUNTERPERALIGNMENT = SEQCOUNTERPERALIGNMENT + 1
				AMINOACIDSPERSTRAIN = AMINOACIDSPERSTRAIN + 1
				if character in 'G':
					G = G + 1
				elif character in 'A':
					A = A + 1
				elif character in 'L':
					L = L + 1
				elif character in 'M':
					M = M + 1
				elif character in 'F':
					F = F + 1
				elif character in 'W':
					W = W + 1
				elif character in 'K':
					K = K + 1
				elif character in 'Q':
					Q = Q + 1
				elif character in 'E':
					E = E + 1
				elif character in 'S':
					S = S + 1
				elif character in 'P':
					P = P + 1
				elif character in 'V':
					V = V + 1
				elif character in 'I':
					I = I + 1
				elif character in 'C':
					C = C + 1
				elif character in 'H':
					H = H + 1
				elif character in 'R':
					R = R + 1
				elif character in 'N':
					N = N + 1
				elif character in 'D':
					D = D + 1
				elif character in 'T':
					T = T + 1
				elif character in 'Y':
					Y = Y + 1
			elif character in '-':
				GAP = GAP + 1
			else:
				NONE = NONE + 1
		if AMINOACIDSPERSTRAIN != (G + A + L + M + F + W + K + Q + E + S + P + V + I + C + Y + H + R + N + D + T):
			print ("\n\n\n THIS IS NOT RIGHT! \n\n\n")
		pG = G / AMINOACIDSPERSTRAIN
		pA = A / AMINOACIDSPERSTRAIN
		pL = L / AMINOACIDSPERSTRAIN
		pM = M / AMINOACIDSPERSTRAIN
		pF = F / AMINOACIDSPERSTRAIN
		pW = W / AMINOACIDSPERSTRAIN
		pK = K / AMINOACIDSPERSTRAIN
		pQ = Q / AMINOACIDSPERSTRAIN
		pE = E / AMINOACIDSPERSTRAIN
		pS = S / AMINOACIDSPERSTRAIN
		pP = P / AMINOACIDSPERSTRAIN
		pV = V / AMINOACIDSPERSTRAIN
		pI = I / AMINOACIDSPERSTRAIN
		pC = C / AMINOACIDSPERSTRAIN
		pH = H / AMINOACIDSPERSTRAIN
		pR = R / AMINOACIDSPERSTRAIN
		pN = N / AMINOACIDSPERSTRAIN
		pD = D / AMINOACIDSPERSTRAIN
		pT = T / AMINOACIDSPERSTRAIN
		pY = Y / AMINOACIDSPERSTRAIN


		idx = record.id

		id = str(idx)
		CURRENTSPECIES = id_to_species[str(idx)]

	#	file.write ('\nG,' + str(pG) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nA,' + str(pA) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nL,' + str(pL) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nM,' + str(pM) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nF,' + str(pF) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nW,' + str(pW) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nK,' + str(pK) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nQ,' + str(pQ) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nE,' + str(pE) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nS,' + str(pS) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nP,' + str(pP) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nV,' + str(pV) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nI,' + str(pI) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nC,' + str(pC) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nH,' + str(pH) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nR,' + str(pR) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nN,' + str(pN) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nD,' + str(pD) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nT,' + str(pT) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))
	#	file.write ('\nY,' + str(pT) + ',' + str(CURRENTSPECIES) + ',' +  str(id) + ',' + str(FILE))

		file.write ('\n' + str(CURRENTSPECIES) +  ',' + str(id) + ',' + str(FILE) + ',' + str(pA) + ',' + str(pC) + ',' + str(pD) + ',' + str(pE) + ',' + str(pF) + ',' +  str(pG) + ',' +  str(pH) + ',' +  str(pI) + ',' +  str(pK) + ',' +  str(pL) + ',' +  str(pM) + ',' +  str(pN) + ',' + str(pP) +  ',' + str(pQ) +  ',' + str(pR) +  ',' + str(pS) + ',' +  str(pT) + ',' + str(pV) + ',' + str(pW) + ',' +  str(pY))

		#print ('\n There are %s amino acids in this strain for this alignment, excluding gaps.' %AMINOACIDSPERSTRAIN)
		AMINOACIDSPERSTRAIN = 0

		G = 0
		A = 0
		L = 0
		M = 0
		F = 0
		W = 0
		K = 0
		Q = 0
		E = 0
		S = 0
		P = 0
		V = 0
		I = 0
		C = 0
		Y = 0
		H = 0
		R = 0
		N = 0
		D = 0
		T = 0

		#print ('\n The local ID for this alignment is %s \n' %id)


		#print 'The problem is the line below'
		#CURRENTSPECIES = [id_to_species('119718979')]
		#print species_to_id['Bacillus_anthracis_str_Ames']



		#print (' \n These amino acids are from the strain, %s\n' %CURRENTSPECIES)

		#print ('\n ----------------------------------------------------------------------------------- \n')
