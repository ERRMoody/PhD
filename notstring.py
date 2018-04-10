def aminoacidchecker(inputfile):
	#this defines a new function, named aminoacidchecker, which takes one input file 'inputfile'



	AAstring = ("GALMFWKQESPVICYHRNDT")
	#This is a string of the 20 aminos acids

	myString = ("GACICNDDDDDDDDTTTTGALMGALMALGMALGMALGMGTTTTTIVVIVYVIVYVIVYVIVYVIVYVIVYNH")
	#this is an example string of amino acids to test the code

	for character in myString:
		print(character)
	#this prints the character for each present in the string of example AAs

	aa_count = {}
	totalcount = 0
	#this creates an empty dictionary

	for character in AAstring:
		#for each character in the string of all amino acids

		aa_count[character] = 0

		#the amino acid count variable 20 amino acids is 0

	for character in myString:
		if character in aa_count:
			aa_count[character] += 1
			totalcount = totalcount + 1

		else:
			exit
	#print aa_count[aa]
		
	for character in AAstring:
		print ("\n --- " + character + " ---")
		print ("There are " + str(aa_count[character]) + " copies of the amino acid " + str(character) + ", out of a total of " + str(totalcount))
		foo = aa_count[character]
		bar = totalcount
		percentageofAA = float(foo)/float(bar)
		percentageofAA = percentageofAA 	
		print("Proportion of " + character + " present is " + str(percentageofAA) + "")



	#mydict[species][character].append(aa_prop)  		


