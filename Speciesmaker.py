import pandas as pd

df = pd.read_csv('fixedlistall.csv')


percaa = df.pivot_table('Percentage of Amino Acid', ['Species', 'Domain','Actual OGT'], 'Amino Acid')


df3 = percaa.groupby(['Species','Domain','Actual OGT'])['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T'].mean()

print (df3)

#df3.to_csv('Moodyspeciesdata.csv')

