import pandas as pd

df2 = pd.read_csv('finalfixedencompassing.csv')

percaa2 = df2[['Species','Strain', 'Domain','Actual OGT','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']]


percaa2.to_csv('Moodystraindata.csv')
