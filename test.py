
import pandas as pd
import csv
import re

#idandspecieslist = pd.read_csv ('nameid.csv')
#specieslist = pd.read_csv ('specieslist.csv')


id = str(386632499)

species_to_id = {}
id_to_species = {}

f = open('nameid.csv')
for line in f:
  fields = re.split(",", line.rstrip())
  species_to_id[fields[0]] = fields[1]
  id_to_species[fields[1]] = fields[0]

print species_to_id['Bacillus_anthracis_str_Ames']
print id_to_species['386632499']