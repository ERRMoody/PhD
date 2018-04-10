from Bio import SeqIO
record_handle = SeqIO.parse("6.phy", "phylip")
for record in record_handle:
	print("%s %i" % (record.id, len(record)))
	for character in str(record.seq):
		if character in '-ARNDBCEQZGHILKMFPSTWYV':
			continue
		else:
			print "Got non-AA character: " + character	