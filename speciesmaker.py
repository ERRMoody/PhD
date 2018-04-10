import pandas as pd

df = pd.read_csv('finalfixedencompassing.csv')


percaa = df[['Species', 'Domain','Actual OGT','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']]


#df3 = percaa.groupby(['Species','Domain','Actual OGT'])['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T'].mean()

#print (df3)

percaa.to_csv('Moodyspeciesdata.csv')
